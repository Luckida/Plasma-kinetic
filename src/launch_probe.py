"""Standard-library-only runtime check, before importing numerical modules.

JSON output is ASCII so Windows console encodings do not corrupt paths.
The script also parses on older Python versions so they can be rejected cleanly.
"""
import json
import os
import platform
import struct
import sys


def check_tk():
    import tkinter
    interpreter = tkinter.Tcl()
    return str(interpreter.eval("info patchlevel"))


def inspect_runtime():
    result = dict(executable=os.path.abspath(sys.executable),
                  version=platform.python_version(),
                  pointer_bits=struct.calcsize("P") * 8,
                  implementation=platform.python_implementation(),
                  tkinter=False, problems=[])
    if sys.version_info[:2] < (3, 11):
        result["problems"].append("Python 3.11 or newer is required")
    if result["pointer_bits"] != 64:
        result["problems"].append("64-bit Python is required")
    if result["implementation"] != "CPython":
        result["problems"].append("CPython is required by the numerical dependencies")
    try:
        result["tcl_version"] = check_tk()
        result["tkinter"] = True
    except Exception as error:
        result["problems"].append("Tkinter/Tcl is unavailable: " + str(error))
    result["ok"] = not result["problems"]
    return result


if __name__ == "__main__":
    report = inspect_runtime()
    print(json.dumps(report, ensure_ascii=True))
    sys.exit(0 if report["ok"] else 1)
