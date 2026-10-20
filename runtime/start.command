#!/bin/bash
set -eu
PATH=/usr/bin:/bin:/usr/sbin:/sbin
COURSE_RESOURCES="$(cd -- "$(dirname -- "$0")/resources" && pwd)"
COURSE_RUNTIME="$COURSE_RESOURCES/work/runtime"
export COURSE_RUNTIME
COURSE_TARGET="$(cat "$COURSE_RESOURCES/runtime/platform.txt")"
COURSE_LANGUAGE="$(cat "$COURSE_RESOURCES/runtime/language.txt")"
say() { if [ "$COURSE_LANGUAGE" = en ]; then printf '%s\n' "$2"; else printf '%s\n' "$1"; fi; }
COURSE_CPU="$(uname -m)"
if [ "$(sysctl -n hw.optional.arm64 2>/dev/null || true)" = "1" ]; then COURSE_CPU="arm64"; fi
if { [ "$COURSE_TARGET" = macos-arm64 ] && [ "$COURSE_CPU" != arm64 ]; } || { [ "$COURSE_TARGET" = macos-intel ] && [ "$COURSE_CPU" != x86_64 ]; }; then
  say '[ARCH] 请下载与 Mac 芯片匹配的课程包。' '[ARCH] Download the package matching your Mac chip.'
  exit 1
fi
if ! mkdir -p "$COURSE_RESOURCES/work"; then
  say '[WRITE] 请将完整课程文件夹移到自己可写的位置。' '[WRITE] Move the complete course folder to a writable location.'
  exit 1
fi
# macOS ships Perl. The lock descriptor survives exec into the preparation shell
# and is released by the OS when that shell and its children exit.
if [ "${1-}" != --prepare-runtime ]; then
  /usr/bin/perl -e '$^F=1024; open(my $lock, ">>", shift) or die "[WRITE] Cannot open initialization lock\n"; $SIG{ALRM}=sub{die "[BUSY] Another window is preparing the runtime\n"}; alarm 300; flock($lock, 2) or die "Cannot lock runtime\n"; alarm 0; exec @ARGV or die "Cannot start preparation\n";' "$COURSE_RESOURCES/work/.runtime-init.lock" /bin/bash "$0" --prepare-runtime
  exec "$COURSE_RUNTIME/bin/python" -E -s "$COURSE_RESOURCES/launch.py" "$@"
fi
CURRENT_LOCATION=""
if [ -f "$COURSE_RUNTIME/.course-root" ]; then CURRENT_LOCATION="$(cat "$COURSE_RUNTIME/.course-root")"; fi
if [ "$CURRENT_LOCATION" != "$COURSE_RESOURCES" ]; then
  say '正在检查并准备内置环境，首次打开需要稍等片刻。' 'Checking and preparing the bundled environment. First launch may take a moment.'
  EXPECTED="$(cat "$COURSE_RESOURCES/runtime/sha256.txt")"
  ACTUAL="$(shasum -a 256 "$COURSE_RESOURCES/runtime/seed.tar.gz" | awk '{print $1}')"
  if [ "$ACTUAL" != "$EXPECTED" ]; then say '[PACKAGE] 环境文件校验失败，请重新解压完整课程包。' '[PACKAGE] Runtime checksum failed. Extract a fresh complete package.'; exit 1; fi
  FREE_KB="$(df -Pk "$COURSE_RESOURCES/work" | awk 'NR==2 {print $4}')"
  NEED_KB="$(cat "$COURSE_RESOURCES/runtime/min_free_kb.txt")"
  if [ "$FREE_KB" -lt "$NEED_KB" ]; then say '[SPACE] 磁盘空间不足，请释放空间后重试。' '[SPACE] Not enough free disk space. Free space on this drive and retry.'; exit 1; fi
  STAGE="$(mktemp -d "$COURSE_RESOURCES/work/runtime-new.XXXXXXXX")"
  # Only generated environments are replaced. Notebook and result folders are separate.
  if ! tar -xzf "$COURSE_RESOURCES/runtime/seed.tar.gz" -C "$STAGE"; then rm -rf -- "$STAGE"; exit 1; fi
  rm -rf -- "$COURSE_RUNTIME"
  mv -- "$STAGE" "$COURSE_RUNTIME"
  "$COURSE_RUNTIME/bin/python" -E -s "$COURSE_RUNTIME/bin/conda-unpack"
  printf '%s\n' "$COURSE_RESOURCES" > "$COURSE_RUNTIME/.course-root.tmp"
  mv -- "$COURSE_RUNTIME/.course-root.tmp" "$COURSE_RUNTIME/.course-root"
fi
exit 0
