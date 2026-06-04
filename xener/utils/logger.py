import logging
import sys
from pathlib import Path
from typing import Optional


class UnbufferedFileHandler(logging.FileHandler):
    """
    Unbuffered file log handler that flushes to disk immediately after each write.
    """
    def __init__(self, filename, mode='a', encoding='utf-8', delay=False):
        super().__init__(filename, mode=mode, encoding=encoding, delay=delay)

    def emit(self, record):
        super().emit(record)
        self.flush()


def setup_logger(
    name: str = "xener",
    log_file: Optional[str] = None,
    level: int = logging.INFO,
    format_string: Optional[str] = None,
    unbuffered: bool = False,
    overwrite: bool = False
) -> logging.Logger:
    """
    Configure and return a project-specific logger.

    Args:
        name: Logger name, defaults to "xener".
        log_file: Log file path; if None, output only to console.
        level: Log level, defaults to logging.INFO.
        format_string: Custom log format string.
        unbuffered: If True, logs are written to disk in real time.
        overwrite: If True, the log file is cleared on each run; defaults to False (append mode).

    Returns:
        Configured logging.Logger object.
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.handlers.clear()

    if format_string is None:
        format_string = "[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s"

    formatter = logging.Formatter(format_string, datefmt="%Y-%m-%d %H:%M:%S")

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_mode = 'w' if overwrite else 'a'
        if unbuffered:
            file_handler = UnbufferedFileHandler(log_file, encoding="utf-8", mode=file_mode)
        else:
            file_handler = logging.FileHandler(log_file, encoding="utf-8", mode=file_mode)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def add_file_handler(
    logger: Optional[logging.Logger] = None,
    log_file: Optional[str] = None,
    level: int = logging.INFO,
    format_string: Optional[str] = None,
    unbuffered: bool = True,
    disable_console: bool = False,
    overwrite: bool = False
) -> logging.Logger:
    """
    Dynamically add a file handler to an existing logger.

    Args:
        logger: Logger object; if None, retrieves the default xener logger.
        log_file: Log file path.
        level: Log level, defaults to logging.INFO.
        format_string: Custom log format string.
        unbuffered: If True, logs are written to disk in real time; defaults to True.
        disable_console: If True, only output to file (console disabled).
        overwrite: If True, the log file is cleared on each run; defaults to False (append mode).

    Returns:
        logging.Logger object with the file handler added.
    """
    if logger is None:
        logger = get_logger()

    if log_file is None:
        raise ValueError("log_file cannot be None")

    if format_string is None:
        format_string = "[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s"

    formatter = logging.Formatter(format_string, datefmt="%Y-%m-%d %H:%M:%S")

    log_path = Path(log_file)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    file_mode = 'w' if overwrite else 'a'
    if unbuffered:
        file_handler = UnbufferedFileHandler(log_file, encoding="utf-8", mode=file_mode)
    else:
        file_handler = logging.FileHandler(log_file, encoding="utf-8", mode=file_mode)
    file_handler.setLevel(level)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    if disable_console:
        remove_console_handler(logger)

    return logger


def remove_console_handler(logger: Optional[logging.Logger] = None) -> logging.Logger:
    """
    Remove the console handler from the logger, keeping only file handlers.

    Args:
        logger: Logger object; if None, retrieves the default xener logger.

    Returns:
        logging.Logger object after removing the console handler.
    """
    if logger is None:
        logger = get_logger()

    console_handlers = [h for h in logger.handlers if isinstance(h, logging.StreamHandler) 
                       and not isinstance(h, logging.FileHandler)]
    for handler in console_handlers:
        logger.removeHandler(handler)

    return logger


def get_logger(name: str = "xener") -> logging.Logger:
    """
    Get a configured logger, creating it with default settings if it does not exist.

    Args:
        name: Logger name.

    Returns:
        logging.Logger object.
    """
    logger = logging.getLogger(name)
    if not logger.handlers:
        return setup_logger(name)
    return logger


class LoggerMixin:
    """
    Logging mixin class that provides convenient logging methods for classes.
    """
    @property
    def logger(self) -> logging.Logger:
        if not hasattr(self, "_logger"):
            self._logger = get_logger(self.__class__.__name__)
        return self._logger

    def debug(self, msg: str, *args, **kwargs):
        self.logger.debug(msg, *args, **kwargs)

    def info(self, msg: str, *args, **kwargs):
        self.logger.info(msg, *args, **kwargs)

    def warning(self, msg: str, *args, **kwargs):
        self.logger.warning(msg, *args, **kwargs)

    def error(self, msg: str, *args, **kwargs):
        self.logger.error(msg, *args, **kwargs)

    def critical(self, msg: str, *args, **kwargs):
        self.logger.critical(msg, *args, **kwargs)


logger = get_logger()
