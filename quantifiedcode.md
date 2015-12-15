[![Code Issues](https://www.quantifiedcode.com/api/v1/project/829d7d4b95d04111b38451d14101545d/badge.svg)](https://www.quantifiedcode.com/app/project/829d7d4b95d04111b38451d14101545d)

[quantifiedcode](https://www.quantifiedcode.com/) is a tool that checks code for best practices and finds bugs. It also  can send a PR to a github project to correct issues and it can run checks on PRs. We are discussing under issue [#583](/MDAnalysis/mdanalysis/issues/583) how to best use *quantifiedcode* with MDAnalysis.

## quantifiedcode: MDAnalysis ##
It is now enabled for MDAnalysis:

https://www.quantifiedcode.com/app/project/gh:MDAnalysis:mdanalysis

We can use quantifiedcode in the same way as one uses the spell checker in a text processor: fix the obvious mistakes and use most of it as suggestions.

## Easy issues for new contributors ##
Most issues are not auto-fixable and in any case, even magic PR don't seem to work at the moment (try clicking on the "Fix me" link that appears for some issues), but from the issue list we could easily generate a few **"easy issues"** for people who want to contribute (with an eye towards GSoC and @dldotson's and @richardjgower's efforts to introduce more people at ASU to contributing to OS projects).

So, if you have a moment to spare, go to the URL above, pick something easy that you think ought to be fixed and raise an issue and *label* it with

* **quantifiedcode** (list of all open [issues tagged "quantifiedcode"](https://github.com/MDAnalysis/mdanalysis/labels/quantifiedcode))
* **difficulty-easy**
* **help wanted**
* ... and any other appropriate labels

As the issue title use the *SAME TITLE* as the QC issue (to avoid duplication).

As example see https://github.com/MDAnalysis/mdanalysis/issues/564 --- assume that someone new to MDAnalysis will read the issue so be explicit and use links to the QC issue.

If the issue has the magic "fix it" link enabled, you can also ask quantifiedcode to submit a PR. If this does not work, look at issue [#583](/MDAnalysis/mdanalysis/issues/583).
