# mean-shift-Kalman-filter-tracker

A MATLAB implementation of a mean-shift tracking algorithm combined it with the Kalman Filter in order to improve its performance. 

Mean-shift localizes the target by finding the coloured histogram of the region that has the best match to the target's histogram. It works reasonably well in most cases but can fail when the target is partially or fully hidden behind an obstacle or by a visually similar object. 
Kalman Filter uses recursive Bayesian estimation to combine the measurement (obtained by the detector) with the prediction (obtained with the dynamic model) and produce a better estimation of the target's location, as seen in the gif below.

![demo](demo.gif)

Red circle represents the variance or uncertainty of the measurement, blue circle the variance of the prediction and the green circle the variance of the combined, final localization. We see that the tracker doesn't fail even when the target is occluded because the uncertainty of the measurement is high - thus the final estimation is influenced more by the prediction.  
**This gif represents the debug view of the tracker. During real case tracking, the circles are not drawn and the video is not playing so slowly, frame-by-frame.*

Testing our solution with several dynamic models (random walk, near constant velocity, near constant acceleration) showed improvement in specific scenarios. The near constant acceleration model proved most successful.
Testing was done with the Tracking Evaluation Toolkit provided in the repository and [VOT sequences](http://www.votchallenge.net/challenges.html).

File porocilo.pdf contains the derivation of the Kalman Filter and an extensive testing report in the Slovenian language.
This project was a personal assignment for the Advanced Topics in Computer Vision course at the University of Ljubljana, Faculty of Computer and Information Science. 
