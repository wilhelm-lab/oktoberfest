# Status

[![pipeline status](https://gitlab.lrz.de/compmass/prosit/app/oktoberfest/badges/develop/pipeline.svg)](https://gitlab.lrz.de/compmass/prosit/app/oktoberfest/commits/develop)

# Example

- `/your/folder/` has some RAW files in it and a `msms.txt` from a MaxQuant search on this RAW files.
- run `make all DATA=/your/folder/`


# CI Stuff
```
python -m venv venv
source venv/bin/activate
# add stuff to pyproject.toml
poetry update
```
If you dont add *oktoberfest* as a *Deploy Key* (GITLAB/REPO/SETTINGS/REPOSITORY) to every repository, the git clone within poetry wont work
