Project #3: Building the NRPZ & interpreting small planets

Project description:
Select an exoplanet with both radial velocity and transit data and download the data (see the instruction in “Downloading_RV_transit_data.docx”). I recommend three systems: HD 209458 b, HD 189733 b, and GJ 436 b. Or you can select whichever exoplanet system that you deem interesting/challenging. 

You will need to:
1, Build the NRPZ using the planet-builder we have been working with
•	Download Fe and Mg abundances from the Hypatia catalog
•	Remove any stars that don’t have both Fe AND Mg abundances
•	Calculate the Fe/Mg ratio for the stars
•	Fit the Fe/Mg distribution to a Gaussian
•	Determine the 3-sigma range of the Fe/Mg distribution
•	Use the lower 3-sigma bound to build a density-mass curve from approximately 0.1-10 Earth masses
•	Use the upper 3-sigma bound to build a density-mass curve from approximately 0.1-10 Earth masses
•	These two density-mass curves will be your NRPZ
2, Use the NEA to plot ‘small planets’ (<2 Earth radii) with mass and radius uncertainties of less than 30% on top of it
3, plot Mercury, Venus, Earth, and Mars
4, Select one planet from EACH of the following lists (so you should have three total). Highlight these planets on your figure.

A few bonus opportunities – you are welcome to do all three, and I encourage you to do so if you have time, but I will only award bonus points for one.
1 calculate the probability that each of the planets from #3 are w/in the NRPZ
2 pick one of the following planets with measured host abundances…
3 compare your NRPZ limits with that of Unterborn+ 2023. What are they different?\
4 Do anything else above and beyond

The questions are:
1, What is the planet mass and the uncertainty?
2, What is the planet radius and the uncertainty?
3, What is the planet density and the uncertainty?
4, How do your measurements of mass, radius, and density compare to other exoplanets with similar masses and radii? You can gather information from NEA and make the comparison plot. 
5, How do you measurements of mass and radius compare with the M-R relation from Chen & Kipping (2016)? 

Goals:
1, Be able to use NEA and other tools to check and download exoplanet data.
2, Strengthen the understanding of measuring planet mass and radius based on observables. 
3, Be able to fit data from exoplanet observations.
4, Understand the measurement uncertainty in observations and how it propagates to measured physical quantities. 
5, Understand how density of a planet is measured and how a single data point is placed in the context of a population. 

Deliverables:
1, A 10-min oral presentation led by one of the team members. We will rotate the leading presenter as the semester goes on. 
	A, The 10-min presentation needs to spare ~2 min for questions.
	B, The presentation should cover motivation, methods, results, and conclusions. 
	C, The presentation should address all questions in the project
D, The presenter should acknowledge the contribution of each team member when appropriate. 

2, A written report led by one of the team members. We will rotate the leading writer as the semester goes on. 
A, The report should be less than 5-page and cover motivation, methods, results, conclusions, and contribution statement for each team member. References do not count towards the 5-page limit.
B, The report should clearly lay out all the assumptions and calculations.
C, The report should address all questions in the project.
D, The report should acknowledge the contribution of each team member when appropriate, e.g., specifying leader and contributors for each section of the report. 

Evaluation and rubric:
Your oral presentation and written report will be evaluated by other teams based the following rubric. Your score will be the average of scores from other teams. 
1, 5 pts for the oral presentation.
	A, Is the presentation clear? (1 pt)
	B, Is the presentation engaging? (1 pt)
C, Is the presentation complete, i.e., covering motivation, methods, results, and conclusions (2 pt)
	D, Does the presentation address all the questions from audience? (1 pt)
2, 5 pts for the written report.
A, Is the report complete, i.e., covering motivation, methods, results, and conclusions (2 pt)
B, Does the repot clearly lay out all the assumptions and calculations? (1 pt)
C, Does the report address all questions in the project? (1 pt)
D, Are the results sensible? (1 pt)

The 10 pt peer evaluation will be converted to 5 pt that goes into the peer review score. In addition, your projects, your questions to presenters and your evaluation of other teams, and self-evaluation will be evaluated by instructors (5 pt). Points are given at a 0.5 pt increment. 

Tips for collaborative research:
1, It is a teamwork, so think carefully how you would use the expertise in your team to accomplish the goals. For example, think what resource you need for the project. By looking at the project description, we need expertise in:
	A, Python including modules such as radvel and pylightcurve
	B, Plotting, or even better graphic design
C, Fitting data with a model, e.g., numerical methods such as least square fitting or Bayesian sampling
D, Error propagation, either analytically or numerically
E, Oral presentation
F, PowerPoint, or other presenting tools
G, Academic writing 
H, Miscellaneous
2, Break the project into small parts and tackle each part based on the expertise that the part requires. 
3, Risk management. Looking at the breakdown of the necessary areas of expertise, it might seem daunting even though you have a team of ~5 people with different backgrounds. Here are my suggestions on how to get meaningful results within a reasonable amount of time.
A, Approach #1, order of magnitude estimation of the RV amplitude and the transit depth, which can be translated to mass and radius through equations given in classes. Find a way of analytically translating the RV and transit measurement uncertainties into mass and radius uncertainties, then the uncertainty for density. 
B, Approach #2, use exiting tools to complete the jobs. I gave an example of the EXOFAST online tool, “EXOFASTRecipeRunningASimpleFit.pdf”. You can read the instruction and complete the tasks. 
C, Approach #3, following the instructions in the python notebook that I provided, “Week_03_Measuring_Mass_Radius_Density.ipynb”, to complete the tasks. This approach may be more complicated than the previous two because it involves multiple python packages and Bayesian inference, but it is generally what is done in research. 
You can take one of the approaches, or multiple approaches, for the purpose of: risk management and cross-validation. Discuss with your teammates, and decide what is the best way given your confidence level and the areas of expertise in your team. 
4, Identifying expertise and specific role is important because it helps to appropriately 5cknowledge the contribution in the written report. 
6, Set explicit deadlines and expectations for each other.
7, Meet regularly at predetermined time.
6, Yell for help.
