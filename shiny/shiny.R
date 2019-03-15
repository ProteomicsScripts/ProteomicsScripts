install.packages('rsconnect')

library(rsconnect)

rsconnect::setAccountInfo(name='lars20070',
                          token='AFD4219F214CFBAEE19AB4170DCC8343',
                          secret='<SECRET>')

rsconnect::deployApp('path/to/your/app')

#https://www.shinyapps.io/admin/#/dashboard