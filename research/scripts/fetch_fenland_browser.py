"""Use a fresh standard browser context for the public PMC download workflow."""
from pathlib import Path
from datetime import datetime,timezone
import asyncio
import hashlib
import json
import zipfile
from playwright.async_api import async_playwright

ROOT=Path(__file__).resolve().parents[2]
URL='https://pmc.ncbi.nlm.nih.gov/articles/instance/9904207/bin/NIHMS1866945-supplement-Supplementary_tables.xlsx'
# Playwright resolves its installed browser; no machine-specific executable path.

async def main():
    record={'requested_url':URL,'transport':'standard_headless_browser_fresh_context','started_utc':datetime.now(timezone.utc).isoformat()}
    destination=ROOT/'research/data/raw/Fenland/Fenland_supplementary_tables.xlsx'
    destination.parent.mkdir(parents=True,exist_ok=True)
    try:
        async with async_playwright() as p:
            browser=await p.chromium.launch(headless=True)
            record['browser_version']=browser.version
            context=await browser.new_context(accept_downloads=True)
            page=await context.new_page()
            async with page.expect_download(timeout=120000) as event:
                try:
                    response=await page.goto(URL,wait_until='domcontentloaded',timeout=45000)
                    if response:record['initial_http_status']=response.status
                    print('Page:',await page.title(),flush=True)
                except Exception as exc:
                    if 'Download is starting' not in str(exc):raise
            download=await event.value
            partial=destination.with_suffix('.xlsx.part')
            await download.save_as(partial)
            with partial.open('rb') as stream:assert stream.read(2)==b'PK','Not an XLSX ZIP'
            with zipfile.ZipFile(partial) as z:
                assert z.testzip() is None
                assert 'xl/workbook.xml' in z.namelist()
            partial.replace(destination)
            record.update(status='downloaded',suggested_filename=download.suggested_filename,
                          path=str(destination.relative_to(ROOT)),bytes=destination.stat().st_size,
                          sha256=hashlib.sha256(destination.read_bytes()).hexdigest())
            await browser.close()
    except Exception as exc:
        record.update(status='not_acquired',error=str(exc))
    record['finished_utc']=datetime.now(timezone.utc).isoformat()
    (ROOT/'research/records/advance_20261017/fenland_browser_retrieval.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2),flush=True)
    return 0 if record['status']=='downloaded' else 1

if __name__=='__main__':raise SystemExit(asyncio.run(main()))
