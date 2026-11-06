# Explicitly run when the environmental library changes; not part of each fit.
args <- commandArgs(FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- Sys.getenv("CAST_EXAMPLE_ROOT", if (length(script)) dirname(script[1L]) else getwd())
source(file.path(root, "inventory.R"), encoding = "UTF-8")
glo_inventory_build(root, Sys.getenv("CAST_DATA_ROOT", "E:/CHECO26/N-SDM/data"))
