Version 0.11.0 of MDAnalysis brought many changes to the API, meaning that many existing scripts will not work.  But panic not, along with these changes, a tool for fixing old scripts, "ten2eleven", has been created.  This page outlines how to use this tool.  For a full list of changes see the the [Release Notes for 0.11.0](ReleaseNotes0110) and for a comprehensive list of changes to be made to code see the [MDAnalysis 0.11 unifying release user guide](MDAnalysis-0.11-unifying-release-user-guide).

### Using ten2eleven

The ten2eleven tool is distributed as part of the MDAnalysis package starting from 0.11.0.
To convert a script, simply call the `ten2eleven` function and pass it the path to the script you want to convert.  This function can accept many filenames.


``` python
import MDAnalysis as mda

mda.ten2eleven('myscript.py', 'myotherscript.py')

```

This will then scan and correct the script in place, after creating a backup called `OLDNAME.bak`.