"""Developer credit shown in the HTML report and the command-line banner.

The logo is not part of the source code: put the official image file in branding/ of the bundle
(bash build_dmdseek.sh --logo FILE) or point DMDSEEK_LOGO to it; it is embedded in every report.
Without it only the text is shown.
"""
import os
import sys
from pathlib import Path

INSTITUTION = 'Seoul National University Hospital'
DEPARTMENT = 'Department of Laboratory Medicine'
LAB = 'Molecular Diagnostics Laboratory'
AFFILIATION = ('Department of Laboratory Medicine, Seoul National University Hospital; '
               'Seoul National University Cancer Research Institute')
PEOPLE = (('Principal investigator', 'Moon-woo Seong, MD, PhD'),
          ('Developer', 'Seung Won Chae'))


def bundle_dir():
    """Folder of the frozen executable, or the source tree (development)."""
    if getattr(sys, 'frozen', False):
        return Path(sys.executable).resolve().parent
    return Path(__file__).resolve().parent.parent


def logo_file():
    env = os.environ.get('DMDSEEK_LOGO')
    if env and Path(env).is_file():
        return Path(env)
    for ext in ('svg', 'png', 'jpg', 'jpeg', 'webp'):
        p = bundle_dir() / 'branding' / f'logo.{ext}'
        if p.is_file():
            return p
    return None
