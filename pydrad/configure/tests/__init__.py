import os


def assert_ignore_blanks(string_1, string_2):
    """
    Compare two multiline strings, ignoring lines that are blank
    """
    string_1_stripped = os.linesep.join([l for l in string_1.splitlines() if l])
    string_2_stripped = os.linesep.join([l for l in string_2.splitlines() if l])
    assert string_1_stripped == string_2_stripped
