import re

def get_leaplog_water_ion_numbers(leap_log_dir):

    with open(leap_log_dir, 'r') as f:
        text = f.read()

    pat = re.compile('> solvatebox[\s\S.]+Added (\d+) residues.')
    matches = re.findall(pat, text)
    water_num = matches[0]

    pat = re.compile('> addIonsRand[\s\S.]+(\d+) Na\+ ion[s]? required to neutralize.')
    matches = re.findall(pat, text)
    Na_num = int(matches[0]) if len(matches) > 0 else 0

    pat = re.compile('> addIonsRand[\s\S.]+(\d+) Cl- ion[s]? required to neutralize.')
    matches = re.findall(pat, text)
    Cl_num = int(matches[0]) if len(matches) > 0 else 0

    return water_num, Na_num, Cl_num
