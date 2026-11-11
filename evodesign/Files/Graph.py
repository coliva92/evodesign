import numpy as np
import numpy.typing as npt



def create_graph_file(contacts_map: npt.NDArray[np.float64],
                      filename: str
                      ) -> None:
    # check that the contact map is a square matrix
    assert contacts_map.shape[0] == contacts_map.shape[1]
    # create the edges
    num_residues = contacts_map.shape[0]
    num_edges = 0
    edges_list = []
    for i in range(num_residues):
        for j in range(num_residues):
            if contacts_map[i][j] == 1:
                num_edges += 1
                edges_list.append([i, j])
    # save the graph info into a text file
    with open(filename, "wt", encoding="utf-8") as txt:
        txt.write(f"{num_residues}\n{num_edges}\n")
        for edge in edges_list:
            txt.write(",".join(edge) + "\n")
    return
