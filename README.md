# STARE2-analysis

This repo contains the data on the event ST 200428A, the FRB from SGR 1935+2154, and the code used to localize the event.

candidate_ovro_200428.fil - filterbank data of ST 200428A from OVRO. 
candidate_gdscc_200428.fil - filterbank data of ST 200428A from GDSCC. 
candidate_delta_200428.fil - filterbank data of ST 200428A from Delta. 

The data are in filterbank format, which consists of a header that begins with HEADER START and ends with HEADER END, followed by unsigned shorts corresponding to the power in a given frequency channel at a given time. The data have been normalized by the passband and channels containing RFI have been replaced with the mean power. 

We note that we use a custom filterbank header that includes the following keys in addition to the standard filterbank header keys (which will cause filterbank readers to fail if this is not accounted for): 
MJD_hour - The hour (in MJD-7hr) that the (typically 5 hr long) observing track started. 
MJD_minute - The minute (in MJD-7hr) that the observing track started. 
MJD_second - The second (in MJD-7hr) that the observing track started. 
MJD_start - The integer day (in MJD-7hr) that the observing track started. 
cand_location - The number of time samples after the start of the observing track that the candidate is located. 
end_sample - The index of the last time sample contained in the filterbank file. 
start_sample - The sample of the first time sample contained in the filterbank file.

localize_event_upload.py - Measures the time delay between two signals in bins by correlating them and fitting the peak of the correlation with a Lorentzian function. The mean of the Lorentzian represents the time delay in bins.

STARE2_localize.ipynb - Given the time delays found, this notebook computes the regions of the sky the signal could have come from.
