#' Get gene related information
#'
#' This function retrieves comprehensive gene information including gene symbols,
#' descriptions, and other annotations for a given set of gene identifiers.
#'
#' @param id A vector of gene identifiers (symbols, Entrez IDs, etc.). Default is NULL.
#' @param org Organism specification. Default is "hs" for human.
#' @param unique Logical, whether to return unique results only. Default is FALSE.
#' @param keepNA Logical, whether to keep NA values in results. Default is TRUE.
#' @param hgVersion Human genome version to use, either "v38" or "v19". Default is c("v38", "v19").
#'
#' @return A data frame containing gene information including symbols, descriptions, and annotations.
#'
#' @examples
#' \dontrun{
#' # Get information for specific genes
#' gene_info <- gene.info(id = c("TP53", "BRCA1", "EGFR"))
#' 
#' # Get unique results only
#' gene_info_unique <- gene.info(id = c("TP53", "BRCA1"), unique = TRUE)
#' }
#'
#' @export
gene.info <- function (id = NULL, org = "hs", unique = FALSE, keepNA = TRUE, 
    hgVersion = c("v38", "v19")) 
{
    replace_greek <- function(id) {
        id <- stringr::str_replace_all(id, "\\u03b1", "alpha")
        id <- stringr::str_replace_all(id, "\\u03b2", "beta")
        id <- stringr::str_replace_all(id, "\\u03b3", "gamma")
        id <- stringr::str_replace_all(id, "\\u03b4", "delta")
        id <- stringr::str_replace_all(id, "\\u03b5", "epsilon")
        id <- stringr::str_replace_all(id, "\\u03bb", "lambda")
        id <- stringr::str_replace_all(id, "\\u03ba", "kappa")
        id <- stringr::str_replace_all(id, "\\u03c3", "sigma")
    }
    replace_back <- function(id) {
        id <- stringr::str_replace_all(id, "alpha", "\\u03b1")
        id <- stringr::str_replace_all(id, "beta", "\\u03b2")
        id <- stringr::str_replace_all(id, "gamma", "\\u03b3")
        id <- stringr::str_replace_all(id, "delta", "\\u03b4")
        id <- stringr::str_replace_all(id, "epsilon", "\\u03b5")
        id <- stringr::str_replace_all(id, "lambda", "\\u03bb")
        id <- stringr::str_replace_all(id, "kappa", "\\u03ba")
        id <- stringr::str_replace_all(id, "sigma", "\\u03c3")
    }
    .initial <- function() {
        pos <- 1
        envir <- as.environment(pos)
        assign(".genekitrEnv", new.env(), envir = envir)
    }
    web.url <- "https://genekitr-china.oss-accelerate.aliyuncs.com"
    ensOrg_name_data <- function() {
        .initial()
        utils::data(list = "ensOrg_name", package = "genekitr", 
            envir = .genekitrEnv)
        get("ensOrg_name", envir = .genekitrEnv)
    }
    mapEnsOrg <- function(organism) {
        if (organism == "hg" | organism == "human" | organism == 
            "hsa" | organism == "hs") 
            organism <- "hsapiens"
        if (organism == "mm" | organism == "mouse") 
            organism <- "mmusculus"
        if (organism == "rn" | organism == "rat") 
            organism <- "rnorvegicus"
        if (organism == "dm" | organism == "fly") 
            organism <- "dmelanogaster"
        if (organism == "dr" | organism == "zebrafish") 
            organism <- "drerio"
        if (organism == "bt" | organism == "bovine") 
            organism <- "btaurus"
        if (organism == "ce" | organism == "worm") 
            organism <- "celegans"
        if (organism == "gg" | organism == "chicken") 
            organism <- "ggallus"
        if (organism == "mmu" | organism == "rhesus") 
            organism <- "mmulatta"
        if (organism == "pt" | organism == "chipm") 
            organism <- "ptroglodytes"
        if (organism == "xenopus") 
            organism <- "xtropicalis"
        ensorg <- ensOrg_name_data()
        rm(ensOrg_name, envir = .genekitrEnv)
        check_all <- apply(ensorg, 2, function(x) organism %in% 
            x)
        if (any(check_all)) {
            org <- ensorg %>% dplyr::filter(eval(parse(text = colnames(.)[check_all])) %in% 
                organism) %>% dplyr::pull(latin_short_name)
        }
        else {
            stop("\nCheck the latin_short_name in `genekitr::ensOrg_name`")
        }
        return(org)
    }
    ensAnno <- function(org, download.method = NULL, hgVersion) {
        if (hgVersion == "v38") {
            org <- mapEnsOrg(tolower(org))
        }
        else {
            org <- mapEnsOrg(tolower(org))
            org <- paste0(org, "_v19")
        }
        data_dir <- tools::R_user_dir("genekitr", which = "data")
        sub_dir <- "/info/gene/"
        data_dir <- paste0(data_dir, sub_dir)
        make_dir(data_dir)
        url <- paste0(web.url, sub_dir, org, "_anno.fst")
        destfile <- paste0(data_dir, org, "_anno.fst")
        web_f_size <- check_web_size(url)
        local_f_size <- file.size(destfile)
        if (is.na(local_f_size)) 
            local_f_size <- 0
        genekitr_download(url, destfile, method = download.method, 
            data_dir, web_f_size, local_f_size)
        dat <- suppressMessages(fst::read.fst(destfile))
        invisible(dat)
    }
    gentype <- function(id, data = NULL, org, hgVersion = "v38") {
        org <- mapEnsOrg(org)
        if (is.null(data)) 
            data <- ensAnno(org, hgVersion = hgVersion)
        if ("symbol" %in% colnames(data)) {
            data_symbol_normal <- data$symbol %>% stringi::stri_remove_empty_na()
            data_symbol_lower <- tolower(data_symbol_normal)
            data_symbol_upper <- toupper(data_symbol_normal)
            n_sym <- sum(id %in% data_symbol_normal | tolower(id) %in% 
                data_symbol_lower | toupper(id) %in% data_symbol_upper)
        }
        else {
            n_sym <- 0L
        }
        if ("ncbi_alias" %in% colnames(data)) {
            data_ncbi_alias_normal <- data$ncbi_alias %>% stringr::str_split("; ") %>% 
                unlist() %>% stringi::stri_remove_empty_na()
            data_ncbi_alias_lower <- tolower(data_ncbi_alias_normal)
            data_ncbi_alias_upper <- toupper(data_ncbi_alias_normal)
            n_ala <- sum(id %in% data_ncbi_alias_normal | tolower(id) %in% 
                data_ncbi_alias_lower | toupper(id) %in% data_ncbi_alias_upper)
        }
        else {
            n_ala <- 0L
        }
        if ("ensembl_alias" %in% colnames(data)) {
            data_ensembl_alias_normal <- data$ensembl_alias %>% 
                stringr::str_split("; ") %>% unlist() %>% stringi::stri_remove_empty_na()
            data_ensembl_alias_lower <- tolower(data_ensembl_alias_normal)
            data_ensembl_alias_upper <- toupper(data_ensembl_alias_normal)
            n_e_ala <- sum(id %in% data_ensembl_alias_normal | 
                tolower(id) %in% data_ensembl_alias_lower | toupper(id) %in% 
                data_ensembl_alias_upper)
        }
        else {
            n_e_ala <- 0L
        }
        if ("ensembl" %in% colnames(data)) {
            data_ensembl <- stringr::str_split(data$ensembl, 
                "; ") %>% unlist() %>% stringi::stri_remove_empty_na()
            n_ens <- sum(id %in% data_ensembl)
        }
        else {
            n_ens <- 0L
        }
        if ("entrezid" %in% colnames(data)) {
            data_entrezid <- data$entrezid %>% stringi::stri_remove_empty_na()
            n_ent <- sum(id %in% data_entrezid)
        }
        else {
            n_ent <- 0L
        }
        if ("uniprot" %in% colnames(data)) {
            data_uniprot <- stringr::str_split(data$uniprot, 
                "; ") %>% unlist() %>% stringi::stri_remove_empty_na()
            n_uni <- sum(id %in% data_uniprot)
        }
        else {
            n_uni <- 0L
        }
        if (sum(n_sym, n_ens, n_ent, n_uni, n_ala, n_e_ala) == 
            0) {
            stop("Wrong organism or input id has no match!")
        }
        else {
            check <- which(c(n_sym, n_ens, n_ent, n_uni, n_ala, 
                n_e_ala) %in% max(n_sym, n_ens, n_ent, n_uni, 
                n_ala, n_e_ala))
            typ <- c("SYMBOL", "ENSEMBL", "ENTREZID", "UNIPROT", 
                "SYMBOL", "SYMBOL")[check] %>% unique()
        }
        return(typ)
    }
    getOrder <- function(org, keytype, download.method = NULL, 
        hgVersion) {
        if (hgVersion == "v38") {
            org <- mapEnsOrg(tolower(org))
        }
        else {
            org <- mapEnsOrg(tolower(org))
            org <- paste0(org, "_v19")
        }
        data_dir <- tools::R_user_dir("genekitr", which = "data")
        sub_dir <- "/info/gene/"
        data_dir <- paste0(data_dir, sub_dir)
        make_dir(data_dir)
        url <- paste0(web.url, sub_dir, org, "_", keytype, "_order.fst")
        destfile <- paste0(data_dir, "/", org, "_", keytype, 
            "_order.fst")
        web_f_size <- check_web_size(url)
        local_f_size <- file.size(destfile)
        if (is.na(local_f_size)) 
            local_f_size <- 0
        genekitr_download(url, destfile, method = download.method, 
            data_dir, web_f_size, local_f_size)
        dat <- suppressMessages(fst::read.fst(destfile))
        invisible(dat)
    }
    make_dir <- function(data_dir) {
        if (!dir.exists(data_dir)) {
            tryCatch({
                dir.create(data_dir, recursive = TRUE)
            }, error = function(e) {
                message(paste0("Seems like you cannot access dir: ", 
                  data_dir, "\nPlease spefify a valid dir to save data..."))
                data_dir <- readline(prompt = "Enter directory: ")
                dir.create(data_dir, recursive = TRUE)
            })
        }
    }
    check_web_size <- function(url) {
        web_f_size <- RCurl::getURL(url, nobody = 1L, header = 1L) %>% 
            strsplit("\r\n") %>% unlist() %>% stringr::str_extract("Content-Length.*[0-9]") %>% 
            stringr::str_remove_all("Content-Length: ") %>% stringi::stri_remove_empty_na() %>% 
            as.numeric()
        return(web_f_size)
    }
    genekitr_download <- function(url, destfile, data_dir, method, 
        web_f_size, local_f_size) {
        if (!file.exists(destfile)) {
            if (is.null(method)) 
                method <- getOption("genekitr.download.method")
            if (!is.null(method) && method != "auto") {
                tryCatch(utils::download.file(url, destfile, 
                  quiet = TRUE, method = method, mode = "wb"), 
                  error = function(e) {
                    message(paste0("Auto download failed...\nPlease download via: ", 
                      url, "\nThen save to: ", data_dir, "\n"))
                  })
            }
            else {
                tryCatch(utils::download.file(url, destfile, 
                  quiet = TRUE, method = "auto", mode = "wb"), 
                  error = function(e) {
                    message(paste0("Auto download failed...\nPlease download manually via: ", 
                      url, "\nThen save to: ", data_dir, "\n"))
                  })
            }
        }
        else if (web_f_size != local_f_size) {
            message("Detected new version data, updating...")
            if (!is.null(method) && method != "auto") {
                tryCatch(utils::download.file(url, destfile, 
                  quiet = TRUE, method = method, mode = "wb"), 
                  error = function(e) {
                    message(paste0("Auto download failed...\nPlease download manually via: ", 
                      url, "\nThen save to: ", data_dir, "\n"))
                  })
            }
            else {
                tryCatch(utils::download.file(url, destfile, 
                  quiet = TRUE, method = "auto", mode = "wb"), 
                  error = function(e) {
                    message(paste0("Auto download failed...\nPlease download manually via: ", 
                      url, "\nThen save to: ", data_dir, "\n"))
                  })
            }
        }
    }
    org <- mapEnsOrg(org)
    hgVersion <- match.arg(hgVersion)
    if (is.null(id)) {
        gene_info <- ensAnno(org, hgVersion = hgVersion)
    }
    else {
        all <- ensAnno(org, hgVersion = hgVersion)
        if (all(id %>% stringr::str_detect(., "ENS"))) 
            id <- stringr::str_split(id, "\\.", simplify = T)[, 
                1]
        id <- replace_greek(id)
        keytype <- gentype(id = id, data = all, org = org, hgVersion = hgVersion) %>% 
            tolower()
        if (keytype == "symbol") {
            input_df <- data.frame(input_id = id)
            id2 <- tolower(id)
            order_dat <- getOrder(org, keytype, hgVersion = hgVersion) %>% 
                dplyr::mutate(`:=`(!!keytype, tolower(.[[keytype]]))) %>% 
                dplyr::filter(eval(parse(text = keytype)) %in% 
                  id2) %>% dplyr::arrange(.[[keytype]])
            input_df <- input_df %>% dplyr::mutate(rep = tolower(input_id)) %>% 
                merge(., order_dat, by.x = "rep", by.y = keytype, 
                  all.x = T)
            gene_info <- all[input_df$rnum, ] %>% dplyr::mutate(input_id = input_df$input_id) %>% 
                dplyr::relocate("input_id", .before = everything()) %>% 
                dplyr::arrange(match(input_id, id))
        }
        else {
            order_dat <- getOrder(org, all_of(keytype), hgVersion = hgVersion) %>% 
                dplyr::filter(eval(parse(text = keytype)) %in% 
                  id) %>% dplyr::mutate(`:=`(!!keytype, factor(.[[keytype]], 
                levels = unique(id)))) %>% dplyr::arrange(.[[keytype]])
            gene_info <- all[order_dat$rnum, ] %>% dplyr::mutate(input_id = order_dat[[keytype]]) %>% 
                dplyr::relocate("input_id", .before = everything()) %>% 
                merge(., as.data.frame(id), by.x = "input_id", 
                  by.y = "id", all.y = T) %>% dplyr::arrange(input_id)
        }
        if (!is.null(id)) {
            tomany_id <- names(table(gene_info$input_id))[table(gene_info$input_id) > 
                1]
            tomany_id <- tomany_id[!tomany_id %in% id[duplicated(id)]]
            if (length(tomany_id) > 0 & length(tomany_id) < 3) {
                message(paste0("Some ID occurs one-to-many match, like \"", 
                  paste0(tomany_id, collapse = ", "), "\"\n"))
            }
            else if (length(tomany_id) > 3) {
                message(paste0("Some ID occurs one-to-many match, like \"", 
                  paste0(tomany_id[1:3], collapse = ", "), "\"...\n"))
            }
        }
        if (unique & length(tomany_id) != 0) {
            sub <- gene_info %>% dplyr::filter(input_id %in% 
                tomany_id)
            other <- gene_info %>% dplyr::filter(!input_id %in% 
                tomany_id)
            if (all(c("entrezid", "chr") %in% colnames(gene_info))) {
                uniq_order <- vapply(tolower(tomany_id), function(x) {
                  res <- c()
                  check <- which(tolower(sub$input_id) == x)
                  n_na <- apply(sub[check, ], 1, function(x) sum(is.na(x)))
                  if (min(n_na) != max(n_na)) {
                    res <- check[which.min(n_na)]
                  }
                  else {
                    res <- NULL
                  }
                  if (keytype == "symbol" && "symbol" %in% colnames(gene_info) && 
                    is.null(res)) {
                    sym <- tolower(sub[check, "symbol"])
                    if (any(sym %in% x)) {
                      res <- check[which(sym %in% x)]
                      if (length(res) > 1) {
                        if ("summary" %in% colnames(sub)) {
                          res <- which(!is.na(sub$summary[res]))
                          if (length(res) > 1) {
                            res <- NULL
                          }
                          else {
                            res <- check[res]
                          }
                        }
                        else {
                          res <- NULL
                        }
                      }
                    }
                    else {
                      res <- NULL
                    }
                  }
                  if (is.null(res) | length(res) == 0) {
                    n_ent <- as.numeric(sub[check, "entrezid"]) %>% 
                      stats::na.omit() %>% as.numeric()
                    if (!max(n_ent) == min(n_ent)) {
                      min_n <- which(n_ent %in% min(n_ent))
                      res <- check[min_n]
                      if (length(min_n) > 1) {
                        if (any(grepl("^[0-9]|X|Y.*$", sub[res, 
                          "chr"]))) {
                          real_chr <- which(grepl("^[0-9]|X|Y.*$", 
                            sub[res, "chr"]))
                          if (length(real_chr) > 1) {
                            res <- res[1]
                          }
                          else {
                            res <- res[real_chr]
                          }
                        }
                        else {
                          res <- res[1]
                        }
                      }
                    }
                    else {
                      res <- check
                      if (any(grepl("^[0-9]|X|Y.*$", sub[res, 
                        "chr"]))) {
                        real_chr <- which(grepl("^[0-9]|X|Y.*$", 
                          sub[res, "chr"]))
                        if (length(real_chr) > 1) {
                          res <- res[real_chr[1]]
                        }
                        else {
                          res <- res[real_chr]
                        }
                      }
                      else {
                        res <- check[1]
                      }
                    }
                  }
                  return(res)
                }) %>% as.numeric()
            }
            else {
                uniq_order <- vapply(tolower(tomany_id), function(x) {
                  res <- c()
                  check <- which(tolower(sub$input_id) == x)
                  n_na <- apply(sub[check, ], 1, function(x) sum(is.na(x)))
                  if (min(n_na) == max(n_na) & keytype != "entrezid") {
                    res <- check[order(as.numeric(tolower(sub$entrezid)[check])) == 
                      1]
                  }
                  else if (min(n_na) == max(n_na) & keytype == 
                    "entrezid") {
                    res <- check[order(as.numeric(tolower(sub$input_id)[check])) == 
                      1]
                  }
                  else if (min(n_na) != max(n_na)) {
                    res <- check[which.min(n_na)]
                  }
                  return(res)
                }) %>% as.numeric()
            }
            row.names(sub) <- NULL
            gene_info <- rbind(other, sub[uniq_order, ])
            gene_info <- gene_info[match(id, gene_info$input_id), 
                ]
        }
        else {
            id <- factor(id, ordered = T, levels = unique(id))
            gene_info$input_id <- factor(gene_info$input_id, 
                ordered = T, levels = unique(id))
            gene_info <- gene_info[order(gene_info$input_id), 
                ]
        }
    }
    if (!is.null(id)) {
        if (keytype %in% c("ensembl", "entrezid", "uniprot")) {
            gene_info <- gene_info %>% dplyr::select(!all_of(keytype))
        }
        else {
            gene_info <- gene_info %>% dplyr::relocate(symbol, 
                .after = input_id)
        }
    }
    if (!keepNA) {
        gene_info <- gene_info %>% filter_at(vars(!input_id), 
            any_vars(!is.na(.)))
    }
    if (!is.null(id)) {
        gene_info$input_id <- replace_back(gene_info$input_id)
    }
    gene_info[] <- lapply(gene_info, as.character)
    rownames(gene_info) <- NULL
    return(gene_info)
}
