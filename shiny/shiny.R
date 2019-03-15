# first shiny app
# follow tutorial at https://www.shinyapps.io/admin/#/dashboard

#install.packages('rsconnect')

library(rsconnect)

#rsconnect::setAccountInfo(name='lars20070', token='AFD4219F214CFBAEE19AB4170DCC8343', secret='BXBlMEalP8WEzCErBjwCKhJ33LF8iUCuX9bDdbU0')

rsconnect::deployApp('path/to/your/app')

#http://shiny.rstudio.com/articles/shinyapps.html
