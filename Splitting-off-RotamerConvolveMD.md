As outlined in the [[MoveToGitHub]], we want to make *applications* separate repository within the MDAnalysis organization. (We should also relabel them simply as *tools*.)

We only have [[RotamerConvolveMD]] in the applications directory at the moment. *RotamerConvolveMD* now lives in its own repository [MDAnalysis/RotamerConvolveMD](https://github.com/MDAnalysis/RotamerConvolveMD).

## Procedure ##
Oliver followed [Detach subdirectory into separate Git repository](http://stackoverflow.com/questions/359424/detach-subdirectory-into-separate-git-repository/17864475#17864475) and the help of `git subtree --help` (working on a fresh clone, just for safety...):
```
git init RotamerConvolveMD
git clone git@github.com:MDAnalysis/mdanalysis.git
cd mdanalysis
git checkout -b develop origin/develop
git subtree split --prefix=applications/RotamerConvolveMD -b split_RotamerConvolveMD
```
Pull the subtree into the new repository:
```
cd ../RotamerConvolveMD
git pull ../mdanalysis split_RotamerConvolveMD
git checkout -b master
git remote add origin git@github.com:MDAnalysis/RotamerConvolveMD.git
git push origin master
```
This created the new repository https://github.com/MDAnalysis/RotamerConvolveMD.

Finally, remove the applications from the mdanalysis directory:
```
cd mdanalysis
git rm -r applications
git commit -m 'removed applications'
git push origin develop
``` 

