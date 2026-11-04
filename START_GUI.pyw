"""Windows: двойной щелчок для запуска Plasma kinetics 2.1 без консоли."""
import os
from pathlib import Path
import sys
import traceback


def error_dialog(message):
    try:
        from tkinter import messagebox
        messagebox.showerror("Plasma kinetics 2.1 — ошибка запуска", message)
    except Exception:
        if os.name == "nt":
            import ctypes
            box = ctypes.windll.user32.MessageBoxW
            box.argtypes = [ctypes.c_void_p, ctypes.c_wchar_p, ctypes.c_wchar_p, ctypes.c_uint]
            box.restype = ctypes.c_int
            box(0, message, "Plasma kinetics 2.1", 0x10)


if __name__ == "__main__":
    try:
        if sys.version_info < (3, 11):
            raise RuntimeError("Нужен Python 3.11 или новее. Установите его с python.org и повторите запуск.")
        sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
        from launcher import main
        main()
    except Exception as error:
        details = traceback.format_exc()
        path = Path(__file__).resolve().parent / "gui_error.log"
        try:
            path.write_text(details, encoding="utf-8")
        except OSError:
            pass
        error_dialog(str(error) + "\n\nПодробности: " + str(path))
