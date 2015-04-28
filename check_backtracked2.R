# Check back tracked results 

# Load the data
outresults_file <- "/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/Results/backtracked_results_1000km_5min.Rdata"
load(outresults_file)
bt_results$Hour <- as.numeric(format(bt_results$timestep, "%H"))
bt_results$Day <- as.numeric(format(bt_results$timestep, "%d"))

# All MCS have IDs < 1502
allmcs <- 1:1502

# Which MCS DO have a Generation point?
genpt <- sort(allmcs[bt_results[bt_results$ID <= 1502 & bt_results$class == "Generation","ID"]])

# Which MCS DON"T have a Generation point?
nogenpts <- allmcs[-bt_results[bt_results$ID <= 1502 & bt_results$class == "Generation","ID"]]

# What happens to the MCS with no generation point?
all_ngp_class <- vector("character")
for (ngp in nogenpts){
    # Get earliest record for this ID, and print class
    first_time <- min(bt_results[bt_results$ID == ngp, "timestep"])
    ft_class <- bt_results[bt_results$ID == ngp & bt_results$timestep == first_time, "class"]
    ft_related <- bt_results[bt_results$ID == ngp & bt_results$timestep == first_time, "related"]
    print(paste("ID:",ngp, " Time:",first_time," Class:",ft_class," Related:",ft_related))
    ft_class <- sub("minor parts from ID [0-9]{1,6}", "minor parts", ft_class)
    all_ngp_class <- append(x = all_ngp_class, values = ft_class)
}
table(all_ngp_class)
# Basically, the majority are small parts of a split in the previous timestep - the largest overlap of the split takes the ID of the previous timestep, but the smaller parts need a new ID.
    
# What is going on with the high ID numbers? 
large_id <- sort(unique(bt_results[bt_results$ID > 1502,"ID"]))
all_lid_class <- vector("character")
end_lid_class <- vector("character")
for (lid in large_id){
    first_time <- min(bt_results[bt_results$ID == lid, "timestep"])
    last_time <- max(bt_results[bt_results$ID == lid, "timestep"])
    ft_class <- bt_results[bt_results$ID == lid & bt_results$timestep == first_time, "class"]
    ft_class <- sub("minor parts from ID [0-9]{1,6}", "minor parts", ft_class)
    lt_class <- bt_results[bt_results$ID == lid & bt_results$timestep == last_time, "class"]
    ft_related <- bt_results[bt_results$ID == lid & bt_results$timestep == first_time, "related"]
    # print(paste("ID:",lid, " Time:",first_time," Class:",ft_class," Related:",ft_related))
    all_lid_class <- append(x = all_lid_class, values = ft_class)
    end_lid_class <- append(x = end_lid_class, values = lt_class)
}
100 * table(all_lid_class) / length(all_lid_class)
# 67% (2285) of the large numbers are POP generations
# 31% (1034) are minor parts of splits in the current timestep
100 * table(end_lid_class) / length(end_lid_class)



