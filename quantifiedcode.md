[![Code Issues (master)](https://www.quantifiedcode.com/api/v1/project/b993444b512c4f9ab23f9f20e6af7479/badge.svg)](https://www.quantifiedcode.com/app/project/b993444b512c4f9ab23f9f20e6af7479)

[quantifiedcode](https://www.quantifiedcode.com/) is a tool that checks code for best practices and finds bugs. It also advertises that it can send PR to a github project that correct issues and it can supposedly run checks on PRs.

It is now enabled for MDAnalysis (but the badge above only points to the _master_ branch whereas we really care about _develop_):

https://www.quantifiedcode.com/app/project/b993444b512c4f9ab23f9f20e6af7479?branch=origin%2Fdevelop

We can use quantifiedcode in the same way as one uses the spell checker in a text processor: fix the obvious mistakes and use most of it as suggestions.

Most issues are not auto-fixable and in any case, even magic PR don't work for me (see other email), but from the issue list we could easily generate a few "easy issues" for people who want to contribute (with an eye towards GSoC and @dldotson's and @richardjgower's efforts to introduce more people at ASU to contributing to OS projects).

So, if you have a moment to spare, go to the URL above, pick something easy that you think ought to be fixed and raise an issue and tag it with

* quantifiedcode
* difficulty-easy
* help wanted

As the issue title use the SAME TITLE as the QC issue (to avoid duplication).

As example see https://github.com/MDAnalysis/mdanalysis/issues/564 --- assume that a newbie will read the issue so be explicit and use links to the QC issue.
