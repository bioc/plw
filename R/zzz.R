.onAttach <- function(libname, pkgname)
#	Add HowToPLW to Windows menu
{
    msg <- sprintf(
        "Package '%s' is deprecated and will be removed from Bioconductor
         version %s", pkgname, "3.12")
    .Deprecated(msg=paste(strwrap(msg, exdent=2), collapse="\n"))
    if( .Platform$OS.type == "windows" && .Platform$GUI == "Rgui" ) winMenuAddItem("Vignettes","plw","HowToPLW()")
}
