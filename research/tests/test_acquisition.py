"""Regression cases from the acquisition audit; no public network requests."""
from contextlib import ExitStack
from hashlib import md5, sha256
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch, Mock
import json
import multiprocessing
import unittest

from research.scripts import download_data as data
from research.scripts import download_gwas as gwas
from research.scripts import prepare_references as refs
from research.scripts.acquisition_state import update_record, read_records, file_lock, check_identity


def concurrent_record_writer(path, barrier, identifier):
    barrier.wait(timeout=15)
    for index in range(4):
        update_record(Path(path), dict(id=f"{identifier}_{index}", status="downloaded"))


def lock_worker(path, active, overlap, barrier):
    barrier.wait(timeout=15)
    with file_lock(Path(path)):
        # Deliberately yield execution while the resource is in use.
        import time
        with active.get_lock():
            active.value += 1
            if active.value != 1:
                overlap.value = 1
        time.sleep(.04)
        with active.get_lock():
            active.value -= 1


class AcquisitionTests(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.manifest = self.root / 'research/records/download_manifest.json'
        self.manifest.parent.mkdir(parents=True)
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)
        for mod in (data, gwas, refs):
            self.stack.enter_context(patch.object(mod, 'ROOT', self.root))
        self.stack.enter_context(patch.object(data, 'MANIFEST', self.manifest))
        self.reference_manifest = self.manifest.with_name('reference_manifest.json')
        self.stack.enter_context(patch.object(refs, 'MANIFEST', self.reference_manifest))
        self.row = dict(id='mvp_gwas', group='gwas', url='https://example.org/v1',
                        path='research/data/raw/test.tsv', role='original', expected_bytes='4')
        self.target = self.root / self.row['path']
        self.target.parent.mkdir(parents=True)
        self.payload = b'GOOD'
        self.old = dict(self.row, status='downloaded', sha256=sha256(self.payload).hexdigest(),
                        bytes=4, retrieved_utc='2026-01-01T00:00:00+00:00')
        metadata = self.root / 'research/data/metadata'
        metadata.mkdir()
        (metadata / 'mvp_md5.txt').write_text(f'{md5(self.payload).hexdigest()} test.tsv\n')

    def seed_cache(self):
        self.target.write_bytes(self.payload)
        update_record(self.manifest, self.old)

    def test_cached_identity_unchanged_is_offline(self):
        self.seed_cache()
        before = self.manifest.read_bytes()
        with patch.object(data.requests, 'get', side_effect=AssertionError('network')):
            result = data.acquire(self.row)
        self.assertEqual(result['sha256'], self.old['sha256'])
        self.assertEqual(self.manifest.read_bytes(), before)

    def test_catalog_changes_reject_without_overwriting_provenance(self):
        self.seed_cache()
        before = self.manifest.read_bytes()
        for key, value in [('url', 'https://example.org/v2'), ('role', 'normalized'),
                           ('expected_bytes', '99'), ('path', 'new.tsv')]:
            with self.subTest(key=key), patch.object(data.requests, 'get', side_effect=AssertionError('network')):
                with self.assertRaisesRegex(ValueError, 'Source identity changed'):
                    data.acquire(dict(self.row, **{key: value}))
                self.assertEqual(self.manifest.read_bytes(), before)
                self.assertEqual(self.target.read_bytes(), self.payload)

    def test_cache_hash_and_size_verified(self):
        self.seed_cache()
        self.target.write_bytes(b'BAD!')
        with self.assertRaisesRegex(ValueError, 'Downloaded input changed'):
            data.acquire(self.row)
        self.seed_cache()
        bad_row = dict(self.row, expected_bytes='99')
        update_record(self.manifest, dict(self.old, expected_bytes='99'))
        with self.assertRaisesRegex(ValueError, 'Cached size mismatch'):
            data.acquire(bad_row)

    def run_gwas(self):
        child = Mock()
        child.wait.return_value = 0
        with patch.object(gwas.subprocess, 'Popen', return_value=child):
            return gwas.acquire_all([self.row])

    def test_transport_success_cannot_hide_missing_file(self):
        self.assertEqual(self.run_gwas(), 1)
        self.assertEqual(read_records(self.manifest)[0]['status'], 'failed')
        self.assertFalse(self.target.exists())

    def test_md5_correct_size_wrong_is_failure(self):
        self.row['expected_bytes'] = '5'
        partial = Path(str(self.target) + '.part')
        partial.write_bytes(self.payload)
        self.assertEqual(self.run_gwas(), 1)
        self.assertFalse(self.target.exists())
        self.assertTrue(partial.exists())
        self.assertIn('Size mismatch', read_records(self.manifest)[0]['error'])

    def test_size_correct_md5_wrong_is_failure(self):
        partial = Path(str(self.target) + '.part')
        partial.write_bytes(b'BAD!')
        self.assertEqual(self.run_gwas(), 1)
        self.assertFalse(self.target.exists())
        self.assertIn('MD5 mismatch', read_records(self.manifest)[0]['error'])

    def test_incomplete_control_file_prevents_promotion(self):
        Path(str(self.target) + '.part').write_bytes(self.payload)
        Path(str(self.target) + '.part.aria2').touch()
        self.assertEqual(self.run_gwas(), 1)
        self.assertFalse(self.target.exists())

    def test_successful_partial_promoted_and_recorded(self):
        Path(str(self.target) + '.part').write_bytes(self.payload)
        self.assertEqual(self.run_gwas(), 0)
        self.assertEqual(self.target.read_bytes(), self.payload)
        self.assertTrue(read_records(self.manifest)[0]['official_md5_verified'])

    def test_existing_verified_file_repairs_stale_failure_offline(self):
        self.target.write_bytes(self.payload)
        update_record(self.manifest, dict(self.row, status='failed', error='old timeout'))
        update_record(self.manifest, dict(id='unrelated', status='downloaded'))
        with patch.object(gwas.subprocess, 'Popen', side_effect=AssertionError('network')):
            self.assertEqual(gwas.acquire_all([self.row]), 0)
        record = next(r for r in read_records(self.manifest) if r['id'] == self.row['id'])
        self.assertEqual(record['status'], 'downloaded')
        self.assertNotIn('retrieved_utc', record)
        self.assertIn('verified_utc', record)
        self.assertEqual(len(read_records(self.manifest)), 2)

    def test_existing_bad_size_is_failure_even_with_matching_md5(self):
        self.target.write_bytes(self.payload)
        self.row['expected_bytes'] = '5'
        with patch.object(gwas.subprocess, 'Popen', side_effect=AssertionError('network')):
            self.assertEqual(gwas.acquire_all([self.row]), 1)
        self.assertEqual(self.target.read_bytes(), self.payload)

    def test_gwas_identity_conflict_keeps_record_and_file(self):
        self.seed_cache()
        before = self.manifest.read_bytes()
        with self.assertRaisesRegex(ValueError, 'Source identity changed'):
            gwas.acquire_all([dict(self.row, role='different')])
        self.assertEqual(self.manifest.read_bytes(), before)

    def test_untracked_partial_cannot_be_resumed(self):
        Path(str(self.target) + '.part').write_bytes(b'OLD')
        with patch.object(data.subprocess, 'run', side_effect=AssertionError('network')):
            with self.assertRaisesRegex(ValueError, 'Untracked partial'):
                data.acquire(self.row, transport='aria2')

    def test_interrupted_download_already_has_source_record(self):
        def interrupt(*args, **kwargs):
            record = read_records(self.manifest)[0]
            self.assertEqual(record['url'], self.row['url'])
            self.assertEqual(record['status'], 'downloading')
            Path(str(self.target) + '.part').write_bytes(b'OLD')
            raise KeyboardInterrupt()
        with patch.object(data.requests, 'get', side_effect=interrupt):
            with self.assertRaises(KeyboardInterrupt):
                data.acquire(self.row)
        with self.assertRaisesRegex(ValueError, 'Source identity changed'):
            data.acquire(dict(self.row, url='https://example.org/v2'))

    def test_source_reservation_and_commit_reject_id_collision(self):
        update_record(self.manifest, dict(self.row, status='downloading'), source=self.row)
        changed = dict(self.row, path='another.tsv')
        with self.assertRaisesRegex(ValueError, 'Source identity changed'):
            update_record(self.manifest, dict(changed, status='downloaded'), source=changed)
        self.assertEqual(read_records(self.manifest)[0]['path'], self.row['path'])

    def test_repeated_query_argument_order_is_significant(self):
        with self.assertRaisesRegex(ValueError, 'Source identity changed'):
            check_identity(dict(url='https://example.org/?id=A&id=B'),
                           dict(url='https://example.org/?id=B&id=A'))

    def test_failed_verification_preserves_historical_evidence_and_recovery(self):
        self.seed_cache()
        self.target.write_bytes(b'BAD!')
        self.assertEqual(gwas.acquire_all([self.row]), 1)
        failed = read_records(self.manifest)[0]
        self.assertEqual(failed['retrieved_utc'], self.old['retrieved_utc'])
        self.assertEqual(failed['previous_verification']['sha256'], self.old['sha256'])
        self.assertNotIn('sha256', failed)
        self.target.write_bytes(self.payload)
        self.assertEqual(gwas.acquire_all([self.row]), 0)
        self.assertEqual(read_records(self.manifest)[0]['retrieved_utc'], self.old['retrieved_utc'])

    def test_reference_request_params_are_part_of_identity(self):
        url = 'https://example.org/geo'
        self.target.write_bytes(self.payload)
        record = dict(path=self.row['path'], url=url+'?acc=GSE1&form=text', sha256=self.old['sha256'])
        update_record(self.reference_manifest, record, key='path')
        with patch.object(refs.requests, 'get', side_effect=AssertionError('network')):
            self.assertEqual(refs.fetch(url, self.target, params={'form':'text', 'acc':'GSE1'}), self.payload)
            with self.assertRaisesRegex(ValueError, 'Reference source changed'):
                refs.fetch(url, self.target, params={'form':'text', 'acc':'GSE2'})

    def test_concurrent_updates_preserve_all_ids(self):
        ctx = multiprocessing.get_context('fork')
        barrier = ctx.Barrier(6)
        jobs = [ctx.Process(target=concurrent_record_writer, args=(str(self.manifest), barrier, i)) for i in range(6)]
        for job in jobs:
            job.start()
        for job in jobs:
            job.join(timeout=20)
            if job.is_alive():
                job.terminate()
                job.join()
            self.assertEqual(job.exitcode, 0)
        self.assertEqual(len(read_records(self.manifest)), 24)

    def test_same_destination_serializes_across_processes(self):
        ctx = multiprocessing.get_context('fork')
        active, overlap, barrier = ctx.Value('i', 0), ctx.Value('i', 0), ctx.Barrier(4)
        jobs = [ctx.Process(target=lock_worker, args=(str(self.target), active, overlap, barrier)) for _ in range(4)]
        for job in jobs:
            job.start()
        for job in jobs:
            job.join(timeout=20)
            if job.is_alive():
                job.terminate()
                job.join()
            self.assertEqual(job.exitcode, 0)
        self.assertEqual(overlap.value, 0)


if __name__ == '__main__':
    unittest.main()
