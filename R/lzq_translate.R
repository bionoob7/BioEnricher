#' @export
lzq_translate <- function (sentence, from = "en", to = "zh", appid = "20231122001888718", 
    key = "5GpDqe8F3pmXfnOkEKGQ") 
{
    seed <- sample(seq_len(1e+05), 1)
    sign <- paste0(appid, sentence, seed, key)
    sign <- openssl::md5(sign)
    res <- NULL
    repeat {
        url <- httr::modify_url("http://api.fanyi.baidu.com/api/trans/vip/translate", 
            query = list(q = sentence, from = from, to = to, 
                appid = appid, salt = seed, sign = sign))
        url <- url(url, encoding = "utf-8")
        res <- jsonlite::fromJSON(url)$trans_result[1, 2]
        if (!is.null(res)) {
            break
        }
    }
    return(res)
}
