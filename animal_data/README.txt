BIC parsed and qced PASS1A metadata. 
Release: 07/01/2019.
Comments: Joslin rats only

We provide the metadata in a single large table called merged_dmaqc_data
  A tab-delimited text file and an RData file are vaillable
In the table we have a row for each vial label (the main sample id). Only ids from samples that were sent to CAS are present.
We also provide the dictionary of the columns in merged_column_dictionary.txt.

Except for viallabel, labelid, pid, bid, barcode, and shipping id the column format is:
"w1.w2.w3", where:
w1.w2 is the column type, which is the form from which it was taken (e.g., acute.test)
w3 is the feature/field name (e.g., distance)

Notable useful columns are:
  viallabel - can be used to map CAS data to the metadata
  labelid,bid,pid - higher level ids (order is pid,bid,labelid,vialid)
  animal.registration.sex
  animal.key.anirandgroup - the experiment group of the rat
  animal.key.timepoint - parsed time point in hours
  animal.key.is_control - whether the rat is an untrained control
  acute.test.distance - the achieved distance during the acute test
  acute.test.howlongshock_seconds - parse shock time (seconds) during the acute test
  acute.test.weight - weight after the acute test
  animal.registration.weight - weight at registration
  columns that start with "calculated.variables": weight and sample time calculated features (e.g., sample freeze time)

  