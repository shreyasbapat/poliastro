import os
from os.path import join, dirname, abspath
from datetime import datetime

import alabaster


# Alabaster theme + mini-extension
html_theme_path = [alabaster.get_path()]
extensions = ["alabaster", "sphinx.ext.intersphinx"]

# Paths relative to invoking conf.py - not this shared file
html_static_path = [join("..", "_shared_static")]
html_theme = "alabaster"
html_favicon = "favicon.ico"
html_theme_options = {
    "logo": "logo_trans.png",
    "logo_name": True,
    "logo_text_align": "center",
    "description": "Astrodynamics in Python",
    "github_user": "poliastro",
    "github_repo": "poliastro",
    "travis_button": True,
    "codecov_button": True,
    "analytics_id": "UA-18486793-1",
    "link": "#3782BE",
    "link_hover": "#3782BE",
    # Wide enough that 80-col code snippets aren't truncated on default font
    # settings (at least for bitprophet's Chrome-on-OSX-Yosemite setup)
    "page_width": "1024px",
}
html_sidebars = {
    "**": ["about.html", "navigation.html", "searchbox.html", "donate.html"]
}

# Regular settings
project = "poliastro"
year = datetime.now().year
copyright = "2013 - %d,  Juan Luis Cano Rodríguez and the poliastro development team" % year
master_doc = "index"
templates_path = ["_templates"]
exclude_trees = ["_build"]
source_suffix = ".rst"
default_role = "obj"
