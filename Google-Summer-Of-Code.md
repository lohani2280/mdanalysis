<img src="https://developers.google.com/open-source/gsoc/images/gsoc2016-sun-373x373.png" title="Google Summer of Code 2016" alt="Google Summer of Code 2016" align="right"/>
MDAnalysis is hosting  **[Google Summer of Code 2016](https://developers.google.com/open-source/gsoc/) students**. In short (more details below — please read the whole page if you are interested!):

* The student will be mentored by at least one experienced MDAnalysis core developer and also receive support from Google. 
* The student will contribute code to MDAnalysis.

Start by reading the [MDAnalysis GSoC blog post](http://www.mdanalysis.org/2016/03/08/gsoc2016/).

The application window closes on **March 25, 2016 at 12:00 (PDT)**.

# Eligibility
Eligibility requirements are described in the [GSoC FAQ](https://developers.google.com/open-source/gsoc/faq): 
* You must be at least 18 years of age.
* You must be a **full or part-time student** at an [accredited university](https://developers.google.com/open-source/gsoc/faq#accredited) (or have been accepted as of April 22, 2016).
* You must be eligible to work in the country you will reside in during the program.
* You have not already participated as a Student in GSoC more than twice.
* You must reside in a country that is not currently embargoed by the United States. See [Program Rules](https://developers.google.com/open-source/gsoc/rules) for more information.

# How to Contact MDAnalysis

You are more than welcome to contact us before submitting your application, we will be happy to advise you on most aspects of the process. Getting in touch first is especially recommended if you are planning to apply to work on an original idea, rather than one of our suggestions.

So if you want to introduce yourself, discuss ideas or your application feel free to write a mail to either the [Discussion](http://groups.google.com/group/mdnalysis-discussion) mailing list or the [Development](http://groups.google.com/group/mdnalysis-devel) mailing list.

# Project Ideas 2016

[Project Ideas](https://github.com/MDAnalysis/mdanalysis/wiki/GSoC-2016-Project-Ideas)


# Who are we?

The MDAnalysis development team is friendly, cooperative and relatively informal. We consist of people from a wide range of backgrounds, including students, PhD-candidates, professors and researchers actively using this software.
We will take your work with respect and appreciate the time that you will spend on MDAnalysis since it will help us in both our own and our community's research projects.

## Available mentors

The following devs have expressed interest in mentoring for GSoC:

- [Max](https://github.com/kain88-de)
- [Manuel](https://github.com/mnmelo)
- [Richard](https://github.com/richardjgowers)

# Am I experienced enough?

The answer is generally: **Yes**. The MDAnalysis Summer of Code team values creativity, intelligence and enthusiasm above specific knowledge of the libraries or algorithms we use. We think that an interested and motivated student who is willing to learn is more valuable than anything else.
The range of available projects should suit people with different backgrounds, both with and without prior experience of computational chemistry.

We also value general software design or development experience above specific library or API knowledge. It is absolutely OK if you have never worked with python before. Just make sure to let us know in your application.

At the same time if you have experience with our tools (such as Python or Cython) or other molecular tools make sure to let us know about that as well.

# Our Expectations from Students

**Communication**

- Write a short report for us once a week
- Commit early and commit often! Push to github so that we can see and review your work.
- Actively work on our project timeline and communicate with us during the community bonding period
- Communicate every working day with your mentor. Just say "Hello" if you like. It can be via email, skype, github comments, etc
- If there is a reason why you can't work or can't contact us on a regular basis please make us aware of this.
- If you don't communicate with us regularly we will fail you.

**Midterm and Final evaluations**

- Set a realistic goal for mid-term. If you fail to meet your own goal we are more likely to fail you in the evaluations
- Have some code merged into our develop branch at the end of the summer to pass the final evaluation
- The last point is a hard requirement. Make sure that your time plan includes it.

# How to write a great Application

Firstly, think about your choice of project carefully, you're going to be doing it for a couple of months, so it's important that you choose something you're going to enjoy. Once you've made your mind up:

1. Make sure you've thought about the project and understand what it entails
2. Don't be afraid to come up with original solutions to the problem
3. Don't be afraid to give us lots of detail about how you would approach the project
4. Contact us early! The earlier you contact us the earlier you will be able to get feedback from us to improve your application

Overall, your application should make us believe that you are capable of completing the project and delivering the functionality to our users. If you aren't sure about anything, get in touch with us, we're happy to advise you.

### Requirements for an Application

We are a taking part in GSoC as under the [Python software foundation](https://www.python.org/psf/) umbrella. This means next to the 
requirements listed here the [requirements of the PSF](https://wiki.python.org/moin/SummerOfCode/2016#How_do_I_Apply.3F) apply.

- Introduce yourself on the [mailing list](https://groups.google.com/forum/#!forum/mdnalysis-devel). Tell us what you plan to work on during 
   the summer or what you have already done with MDAnalysis
- Checkout the source and run the tests from the source code (If you don't manage to setup the virtual environment don't be embarrassed to ask for help)
- Fix an [easy bug](https://github.com/MDAnalysis/mdanalysis/issues?q=is%3Aopen+is%3Aissue+label%3ADifficulty-easy). Alternatively you can also write a new tests or update update our documentation.

### Dividing your project

We require that all of our students have at least one commit in our `develop` branch before the end of the summer. The best way to achieve this is to divide your project into small self contained subprojects and plan to merge at least one of them around midterm. 

A good example for this is the project to port MDAnalysis to python 3. Here each extension and submodule give a natural sub project that you can tackle one after the other. For other it might not be that easy to find suitable sub-tasks but you can always ask us on the mailing list for help.

During your summer you'll encounter bugs in MDAnalysis or find code that can be refactored to help you implement your ideas. You can also immediately fix them in the develop branch and help us all out. This has several advantages. All your pull request will only concentrate on specific features and are much better to review. And you'll also get direct feedback from other devs and user during the summer.

Since this is a hard requirement we as mentors will also have an eye on that and check if your proposal incorporates it and also warn you ahead of time during the summer if we see that you might not make it. Communicating with us on a regular basis is vital for that, though. 

### Your application should include answers to the following questions.

- Why are you interested in working with us?
- Have you used MDAnalysis for your research already?
- Do you have any experience programming?
- Do you have any exams during GSoC or plan a vacation during the summer?

## How to estimate time needed for development

To get a feeling for the code and get some experience with our code you can go and tackle some of our easy bugs. Look at the code that you want to change, check if it follows our coding guidelines. Do some research on the API's you want to use, plan what classes you will add and how their public API will look. Write down your algorithms in pseudo code. The better your research is and the better you plan ahead the easier it will be to judge how long a given task will take. For your time estimates you should also consider that you can do less stuff during exams and try to be a bit conservative. If you have never done anything like GSoC before you will tend to underestimate the time to complete a task. We know that giving these estimates is not easy and that also professionals have problems with it. Having a good plan, knowing its weak and strong points will help a lot. 


