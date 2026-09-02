from xmlload_helper import pyarts, recurse
import os

path = pyarts.arts.globals.data.arts_source_dir
if len(path):
    print(f"arts source found at {path} - starting test run")  
    subs = ["tests", "examples", "doc", "python", 'src', 'testdata']
    paths = [os.path.join(path, sub) for sub in subs]
    for p in paths:
        recurse(p)
    print("All XML files read successfully!")
else:
    print("arts source dir not found - skipping test run")
