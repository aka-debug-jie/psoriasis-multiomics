"""Patch Windows I/O and memory reporting only; statistical sources stay unchanged."""
from pathlib import Path
import sys
root=Path(sys.argv[1]);h=root/'src/mmap_read.hpp';s=h.read_text()
s=s.replace('#if defined(__linux__) || defined(__APPLE__)','#if defined(_WIN32)\n  #ifndef NOMINMAX\n  #define NOMINMAX\n  #endif\n  #include <windows.h>\n#elif defined(__linux__) || defined(__APPLE__)')
s=s.replace('    munmap(addr_, size_);','''#if defined(_WIN32)
    if(addr_) UnmapViewOfFile(addr_);
#else
    munmap(addr_, size_);
#endif''');h.write_text(s)
f=root/'src/mmap_read.cpp';s=f.read_text();start=s.index('MappedFile mmap_file(')
s=s[:start]+'''#if defined(_WIN32)
#include <filesystem>
MappedFile mmap_file(const char* filename) {
  MappedFile result{std::string(filename)};
  auto wide = std::filesystem::u8path(filename).wstring();
  HANDLE file = CreateFileW(wide.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr,
                            OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
  if(file == INVALID_HANDLE_VALUE) return result;
  LARGE_INTEGER size;
  if(!GetFileSizeEx(file, &size) || size.QuadPart <= 0 || size.QuadPart % sizeof(float) != 0) {
    CloseHandle(file); return result;
  }
  HANDLE mapping = CreateFileMappingW(file, nullptr, PAGE_READONLY, 0, 0, nullptr);
  CloseHandle(file);
  if(!mapping) return result;
  void* address = MapViewOfFile(mapping, FILE_MAP_READ, 0, 0, 0);
  CloseHandle(mapping);
  if(address) { result.addr_ = address; result.size_ = static_cast<size_t>(size.QuadPart); }
  return result;
}
#else
'''+s[start:]+'\n#endif\n';f.write_text(s)
f=root/'src/perf_timer.cpp';s=f.read_text().replace('#include <sys/resource.h>','''#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <psapi.h>
#else
#include <sys/resource.h>
#endif''');s=s.replace('  struct rusage usage;\n  getrusage(RUSAGE_SELF, &usage);\n  return usage.ru_maxrss * 1024;', '''#if defined(_WIN32)
  PROCESS_MEMORY_COUNTERS counters{};
  if(!GetProcessMemoryInfo(GetCurrentProcess(), &counters, sizeof(counters))) return 0;
  return counters.PeakWorkingSetSize;
#else
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  return usage.ru_maxrss * 1024;
#endif''');f.write_text(s)
f=root/'CMakeLists.txt';s=f.read_text();s+='\nif(WIN32)\n target_link_libraries(smr psapi ws2_32)\nendif()\n';f.write_text(s)
for relative in ['src/mmap_read.hpp','src/mmap_read.cpp','src/perf_timer.cpp','CMakeLists.txt']:
 f=root/relative
 prefix='#'if f.suffix=='.txt'else'//'
 f.write_text(prefix+' Modified 2026-10-21 by the psoriasis-multiomics course: Windows Unicode I/O, memory reporting and system-library linkage.\n'+f.read_text())
