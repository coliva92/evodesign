from argparse import ArgumentParser
import evodesign.Utils as Utils





if __name__ == "__main__":
    parser = ArgumentParser(prog="average_fitness",
                            description="A utility script that computes the average "
                                        "fitness for each population in the given "
                                        "workspace.")
    parser.add_argument("workspace_dir",
                        help="path to the workspace folder for which populations the "
                             "average fitness will be computed")
    parser.add_argument("output_path",
                        help="path to the file where the computed averages will be stored")
    parser.add_argument("-c", "--fitness_column",
                        type=str, 
                        default=None,
                        help="the name of the column in the population CSV files that "
                             "stores the fitness values")
    args = parser.parse_args()
    averages = Utils.workspace_average_fitness(args.workspace_dir, args.fitness_column)
    averages.to_csv(args.output_path)
