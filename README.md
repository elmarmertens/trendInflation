# trendInflation

code for trend inflation model of Mertens (2016, REStat, http://dx.doi.org/10.1162/REST_a_00549) and Garnier, Mertens and Nelson (2015, IJCB, http://www.ijcb.org/journal/ijcb15q4a2.htm).

The model implemented here is the variant with SV in all inflation gaps.

The `master` branch works with US as well as EuroAre data. US data provided by FRED (for realized and trimmed inflation; use of surveys is pending). Euro area data as obtained from the ECB's SDW, including Euroarea SPF expectations.

Complete replication files for the paper are here: http://dx.doi.org/10.7910/DVN/MUF2HC

## TODO
- construct data sets with US-SPF
- add arguments: Nsim, burnin, p etc.
- add SV(O)-(t) shocks
- Minnesota priors, CTA or equationfilter code for VAR?
- consider SV ratio process?
