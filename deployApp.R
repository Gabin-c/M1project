library(BiocManager)
options(repos = BiocManager::repositories())

rsconnect::setAccountInfo(name='gabin-coudray', token='86AE0832EDD4932875AED002854FC229', secret='EL3cI7L3uQo0Y3TO8djCZn1BBAyT+twk/LQU+WeU')
library(rsconnect)
rsconnect::deployApp('/Users/Gabin/Desktop/M1Project/M1Project/shiny/M1Project')


