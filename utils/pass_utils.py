

def remove_invalid_characters(s):
    return s.replace(".", "")

def parse_pass_spectra(s, ):
    split = s.split()
    pa = split[0]
    pi = split[1]
    activity = " ".join(split[2:])
    return activity, {"Pa": pa, "Pi": pi}