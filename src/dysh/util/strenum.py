# thanks Claude!
from enum import Enum


class StrEnum(str, Enum):
    """
    Enum where members are also (and must be) strings.

    This is a backport of Python 3.11's StrEnum for use in Python 3.10.
    """

    def __new__(cls, value):
        if not isinstance(value, str):
            raise TypeError(f"{value!r} is not a string")
        member = str.__new__(cls, value)
        member._value_ = value
        return member

    def __str__(self):
        return self._value_

    def _generate_next_value_(name, start, count, last_values):
        """
        Return the lower-cased version of the member name.
        """
        return name.lower()
