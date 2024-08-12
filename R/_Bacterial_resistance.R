# -------------------------------------------------------------------------
# A classification system for bacterial resistance
# Cesar A. Saavedra
# Estadistico, Universidad del Valle
# Cali, Colombia
# 2024
# -------------------------------------------------------------------------
## Load packages
suppressMessages(library(openxlsx))
suppressMessages(library(readxl))
suppressMessages(library(data.table))
suppressMessages(library(openxlsx))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(lubridate))
# -------------------------------------------------------------------------
## R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
root <- "/Users/cesara.saavedravanegas/Documents/Investigacion_Javeriana" # Ruta de acceso directorio de trabajo
# -------------------------------------------------------------------------
## Datos Whonet

# -------------------------------------------------------------------------
## Families
# Aeromonadales
# Bacillales
# Enterobacteriales
# Lactobacillales
# Pseudomonadales
# Saccharomycetales
# Staphylococcus epidermidis

## Levels
# Intermediate = I
# Sensitive = S
# Resistant =R

# -------------------------------------------------------------------------
## Split by family to asociate category of resistance 


# -------------------------------------------------------------------------