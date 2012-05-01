#' Build raster filename
#' 
#' Returns the proper filename to use with raster given the base filename and the file format.
#' 
#' @param filename_base Character. The "output filename" used by writeRaster (the base filename).
#' @param format Character. The output file type.  Default is "raster".
#' @return Returns a character string that should be usable by raster(), brick(), or stack()
#' @author Jonathan A. Greenberg and Robert Hijimans
#' @seealso \code{\link[raster]{writeFormats}}, \code{\link[raster]{writeRaster}}
# @keywords {weighted covariance matrix}
# {weighted means}
# @examples
# \dontrun{
# } 
#' @export

build_raster_filename <- function(filename_base,format="raster")
{
	full_filetype = raster:::.filetype(format=format, filename=filename_base)
	full_filename = raster:::.getExtension(filename_base, full_filetype)
	return(full_filename)	
}