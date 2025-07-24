#' @export
lzq_quiet <- function (..., messages = FALSE, cat = FALSE) 
{
    if (!cat) {
        sink(tempfile())
        on.exit(sink())
    }
    out <- if (messages) {
        eval(...)
    }
    else {
        suppressMessages(eval(...))
    }
    out
}
