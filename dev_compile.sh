#!/bin/sh
cp setup.py setup_tmp.py
cp dev_setup.py setup.py
pip install . -v
rm setup.py
mv setup_tmp.py setup.py