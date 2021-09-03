library(ggplot2)
options(bitmapType = "cairo-png")
pheno_df <- read_tsv("data/all_vars.tsv") %>%
    filter(pop != "EUR") %>%
    mutate(sex = c("Female", "Male")[sex + 1])

icd10_names <- trait_df %>%
    filter(phenotype %in% icd10_codes) %>%
    select(code = phenotype, description) %>%
    mutate(description = substr(description, 29, 100))

### Exploratory plots

p1 <- ggplot(pheno_df, aes(pop)) +
    geom_bar(aes(fill = sex, group = sex))

ggsave("plots/count.png", p1, height = 5)

icd10_plot_df <- pheno_df %>%
    select(pop, sex, icd10_names$code) %>%
    group_by(pop, sex) %>%
    summarise(across(everything(), ~ sum(.x, na.rm = T))) %>%
    tidyr::pivot_longer(-c(pop, sex), names_to = "code", values_to = "n") %>%
    left_join(icd10_names, by = "code")

p2 <- ggplot(icd10_plot_df, aes(pop, n)) +
    geom_bar(stat = "identity", aes(fill = sex, group = sex)) +
    facet_wrap("description", scales = "free", labeller = label_wrap_gen())

ggsave("plots/icd10.png", p2, width = 10)

noncancer
non_icd10_id <- c(
    "dvt" = "6152_5", "nopain" = "6159_100", "NC1111" = "20002_1111",
    "NC1226" = "20002_1226", "NC1473" = "20002_1473", "NC1466" = "20002_1466",
    "NC1265" = "20002_1265", "NC1309" = "20002_1309", "NC1452" = "20002_1452",
    "NC1162" = "20002_1162"
)

non_icd10_names <- trait_df %>%
    filter(phenotype %in% c(non_cancer_codes, "6152_5", "6159_100")) %>%
    select(phenotype, description)

non_icd10_plot_df <- pheno_df %>%
    select(pop, sex, names(non_icd10_id)) %>%
    group_by(pop, sex) %>%
    summarise(across(everything(), ~ sum(.x, na.rm = T))) %>%
    tidyr::pivot_longer(-c(pop, sex), names_to = "code", values_to = "n") %>%
    mutate(phenotype = non_icd10_id[code]) %>%
    left_join(non_icd10_names, by = "phenotype")

p3 <- ggplot(non_icd10_plot_df, aes(pop, n)) +
    geom_bar(stat = "identity", aes(fill = sex, group = sex)) +
    facet_wrap("description", scales = "free", labeller = label_wrap_gen())

ggsave("plots/non_icd10.png", p3, width = 10)