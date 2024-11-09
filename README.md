# Wide-Field-Inference-2024
Inferring spiking rate from wide-field recorded fluorescence.

This repository implements the inference solution described in the following: 
https://www.biorxiv.org/content/10.1101/2024.10.18.619111v1.abstract

To start, follow one of the examples according to the scenario that fits your data best:
A. You have long (hour) recordings and 
(i) You wish to maintain as much information as possible in the inferred spiking rates. infer_neural_from_long_dffed_noisy_fluorescence.m
(ii) You wish to smooth the results. another_example_fluorescence_audio_gcamp6.m
B. You have multiple short trial recordings and 
(i) You wish to maintain as much information as possible in the inferred spiking rates. infer_neural_from_short_dffed_noisy_fluorescence.m
(ii) You wish to smooth the results*. (* no example is available yet for this case)

To appropriately run one of these examples using your data, you should know gamma - the amount your calcium indicator decays in a single bin.
For example, in 10hz recordings with GCaMP6s gamma = 0.95;
In 40hz recordings with GCaMP6f gamma = 0.97;

You are welcome to send data and Inquiries. 

Unpublished data is courtesy of the McCormick lab at the University of Oregon and the Zagha lab at the University of California Riverside.

This site is under construction. 
