import networkx as nx
import pandas as pd
import collections
import argparse as ap
import os


def parse_gfa(file_path):
    """Parses a GFA file and returns a graph."""
    G = nx.Graph()
    
    with open(file_path, 'r') as f:
        for line in f:
            
            parts = line.strip().split("\t")

            if parts[0] == "S":  # Segment line
                G.add_node(parts[1])
            elif parts[0] == "L":  # Link line
                G.add_edge(parts[1], parts[3])  # From -> To

    return G

def main():
    global args
    parser = ap.ArgumentParser()
    parser.add_argument("-i", "--gfa", dest="gfa", required=True, help='What is the path to the input gfa file')
    parser.add_argument("-o", "--output", dest='output', required=True, help='Give the full path to output the dataframe')
    parser.add_argument("-s", "--sample", dest='sample', required=False, help='Sample name')
    parser.add_argument("-m", "--mashmap", dest='mashmap', required=False, help='Mashmap file')
    parser.add_argument("-r", "--rdna", dest='rdna', required=False, help='rDNA nodes file')
    parser.add_argument("-p", "--paths", dest='paths', required=False, help='Paths file')
    parser.add_argument("-f", "--scfmap", dest='scfmap', required=False, help='scfmap file')
    parser.add_argument("-c", "--scaffold", dest='scaffold', required=False, help='Scaffold file')
    args = parser.parse_args()

    sample = str(args.gfa).split("/")[-1].split(".")[0]
    if args.sample:
        sample = args.sample
    
    # mySCFMapFile = '/projects/ch-lee-lab/HGSVC/HPRC_Release2/GraphFiles/'+sample+".assembly.scfmap"
    # myMashMapFile = '/projects/ch-lee-lab/HGSVC/HPRC_Release2/GraphFiles/'+sample+".assembly-ref.comp.mashmap"
    # myPathsFile = '/projects/ch-lee-lab/HGSVC/HPRC_Release2/GraphFiles/'+sample+".assembly.paths.tsv"
    input_dir  = "/".join(args.gfa.split("/")[:-2])
    scfMapFile = input_dir + "/" + sample + ".scfmap"
    pathsFile  = input_dir + "/" + sample + ".paths.tsv"
    if not os.path.exists(scfMapFile):
        scfMapFile = input_dir + "/assembly.scfmap"
        pathsFile  = input_dir + "/assembly.paths.tsv"
    
    if not os.path.exists(scfMapFile):
        input_dir  = "/".join(args.gfa.split("/")[:-1])
        scfMapFile = input_dir + "/assembly.scfmap"
    
    if args.paths:
        pathsFile = args.paths
    
    input_dir = "/".join(args.gfa.split("/")[:-1])
    mashMapFile   = input_dir + "/" + sample + ".ref.comp.mashmap"
    rDNAnodesFile = input_dir + "/" + sample + ".rdna.nodes"
    if not os.path.exists(mashMapFile):
        mashMapFile   = input_dir + "/assembly-ref.comp.mashmap"
        rDNAnodesFile = input_dir + "/rdna.nodes"

    if args.mashmap:
        mashMapFile = args.mashmap
    
    if args.rdna:
        rDNAnodesFile = args.rdna
    
    if args.scfmap:
        scfMapFile = args.scfmap
    
    if args.scaffold:
        scaffFile = args.scaffold
    else:
        scaffFile = "T2TScaffolds3.csv"

    # goodArangDF = pd.read_csv("/projects/ch-lee-lab/HGSVC/HPRC_Release2/scripts/T2TScaffolds.csv")
    YscaffDF = pd.read_csv(scaffFile, sep='\t') # Confidently assigned chrY scaffold
    YscaffDict={x:[] for x in YscaffDF['Genome']}
    for genome, seq in zip(YscaffDF['Genome'], YscaffDF['Seq']):
        if seq[:5] == 'chrY_':
            seq = seq.split("chrY_")[1]
        YscaffDict[genome].append(seq)
        
    ## Read in the scf/Paths File
    scfDF = pd.read_csv(scfMapFile, header=None)
    scfScaffoldToPathDict={}
    scfPathToScaffoldDict={}
    for row in scfDF[0]:
        col = row.split(" ")
        if 'path' == col[0]:
            scfScaffoldToPathDict[col[1]] = col[2]
            scfPathToScaffoldDict[col[2]] = col[1]

    pathDF = pd.read_csv(pathsFile, sep='\t')
    # print(pathDF.head())

    ## Collect nodes that are on the Yscaff paths
    YscaffPaths = []
    YscaffPathNodes=[]
    YscaffHap = "hap"
    for scaffold in YscaffDict[sample]:
        if scaffold not in scfScaffoldToPathDict:
            print(f"Skipping scaffold {scaffold} as it is not in the scfScaffoldToPathDict")
            continue
        myPath = scfScaffoldToPathDict[scaffold]
        nodes=list(pathDF[pathDF['name'] == myPath]['path'])[0]
        print(f"Yscaff {scaffold} : {myPath} {nodes}")
        YscaffHap = myPath[:3]
        YscaffPaths.append(myPath)
        for tig in nodes.split(","):
            if tig[:4] == 'utig':
               YscaffPathNodes.append(tig[:-1])
            else:
                continue
    
    print(f"YscaffHap: {YscaffHap}")
    print("YscaffPathNodes: " + ",".join(YscaffPathNodes))
                
    #Read in the mashmap and assign a list of chromosome IDs to nodes
    mashDict={}
    mashDF=pd.read_csv(mashMapFile, sep='\t', header=None)
    for x,y in zip(mashDF[0], mashDF[5]):
        if x in mashDict:
            mashDict[x].append(y)
        else:
            mashDict[x]=[y]
    
    # Read in the rDNA nodes
    rDNA_nodes = []
    with open(rDNAnodesFile, 'r') as f:
        for line in f:
            rDNA_nodes.append(line.strip())
    
    #Which nodes are only chrY or mostly identified as chrY
    mashmapYNodes=[]    # nodes not in YscaffPathNodes but with best hit to chrY via MashMap - check for rDNA components later
    utigToChrDict={}    # final assignment of chromosome IDs to nodes - not used anymore
    xChromosomeNodes=[] # nodes with best hit to chrX via MashMap
    for x,y in mashDict.items():
        tempCount = collections.Counter(y)
        maxID = max(tempCount, key=tempCount.get)
    
        #Look for the nodes on the paths in the paths.tsv first
        if x in YscaffPathNodes:
            utigToChrDict[x]='chrY'
        
        elif x in rDNA_nodes:
            utigToChrDict[x]='rDNA'
        
        else:
            utigToChrDict[x]=maxID
        
            if maxID == 'chrY':
                mashmapYNodes.append(x)
            
            elif maxID == 'chrX':
                xChromosomeNodes.append(x)

            
    # Read graph from GFA file
    gfa_file = str(args.gfa)
    graph = parse_gfa(gfa_file)

    # Iterate through the nodes connected to the YscaffPathNodes
    # and collect all unplaced utig nodes not in YscaffPathNodes nor xChromosomeNodes
    
    visited = set()  # To keep track of visited nodes
    unplacedTigs=[]  # To keep track of unplaced utigs
    for node in YscaffPathNodes:

        if node in visited:
            continue
        
        stack = [node]  # Use a stack to manage the nodes to visit

        while stack:
            current_node = stack.pop()  # Get the last node added to the stack

            if current_node not in visited:
                # print(f"Visiting tig: {current_node}")
                if current_node in xChromosomeNodes:
                    print(f"Skipping tig {current_node} as it is in xChromosomeNodes")
                elif current_node in utigToChrDict and utigToChrDict[current_node] != 'chrY':
                    print(f"Skipping tig {current_node} as it is not chrY")
                else:
                    if current_node not in YscaffPathNodes:
                        # print(f"Adding tig {current_node} to unplaced")
                        # add node to unplacedTigs if neighbor is not in xChromosomeNodes
                        # this is to prevent adding nodes that  belong to PAR1
                        if all(neighbor not in xChromosomeNodes for neighbor in graph.neighbors(current_node)):
                            unplacedTigs.append(current_node)

                    for neighbor in graph.neighbors(current_node):
                        if neighbor[:4] == 'utig' and neighbor not in visited:
                            stack.append(neighbor)

                visited.add(current_node)
    
    if len(unplacedTigs) == 0:
        print("No unplaced tigs found")
    else:
        print(sample + " unplacedTigs: " + ",".join(unplacedTigs))
    # Look for unplaced tigs with hits to Y from mashmapYNodes, exclude if it is within a rDNA component
    
    visited = set() # clear off the visited nodes to traverse the graph
    unplacedTigCandidates = set()  # temporal list to hold Y candidates while traversing
    mashmapUnplaced = []    # list of unplaced tigs with mashmap hits to chrY, not connected to the 'confident Y'
    nodesToSkip = set()
    for node in mashmapYNodes:
        if node in mashmapUnplaced or node in nodesToSkip or node in unplacedTigs:
            continue
        else:
            # Traverse the graph to find all nodes in the component
            stack = [node]
            while stack:
                current_node = stack.pop()
                if current_node not in visited:
                    if current_node in utigToChrDict:
                        if utigToChrDict[current_node] == 'chrY':
                            unplacedTigCandidates.add(current_node)
                        else:
                            nodesToSkip.update(unplacedTigCandidates)
                            unplacedTigCandidates.clear()
                            break # break the loop for traversing the graph
                    else:
                        # Node has no good mashmap hit - possibly a small or repetitive node
                        if current_node == 'rDNA':
                            nodesToSkip.add(current_node)
                            unplacedTigCandidates.clear()
                            break
                        else:
                            unplacedTigCandidates.add(current_node)

                    for neighbor in graph.neighbors(current_node):
                        if neighbor[:4] == 'utig' and neighbor not in visited:
                            stack.append(neighbor)
                visited.add(current_node)
            
            if len(unplacedTigCandidates) > 0:
                mashmapUnplaced.extend(unplacedTigCandidates)
                unplacedTigCandidates.clear()
    
    if len(mashmapUnplaced) == 0:
        print("No unplaced component tigs found")
    else:
        print(sample + " mashmapUnplaced: " + ",".join(mashmapUnplaced))


    # For each of our 'unplaced tigs', see if there is a path that uses it and print the scaffold name
    newList=[]
    isUnplaced  = False
    hasUnplaced = False
    isComponent = False
    listUnplaced = []
    for row in pathDF.itertuples():
        isUnplaced = True
        hasUnplaced = False
        isComponent = False
        listUnplaced.clear()
        if row.name in YscaffPaths:
            newList.append([sample,scfPathToScaffoldDict[row.name],'SCAFFOLD',row.name, row.path])
        else:
            for node in row.path.split(","):
                node=node[:-1]
                if node[:4] != 'utig':
                    continue
                
                # print(f"Checking node {node}")

                if node in unplacedTigs:
                    print(f"Node {node} is in unplacedTigs - {row.name}")
                    listUnplaced.append(node)
                    if row.name in scfPathToScaffoldDict.keys():
                        if row.name[:3] == 'mat' and YscaffHap == 'pat' or row.name[:3] == 'pat' and YscaffHap == 'mat':
                            print(f"Skipping {row.name} as it is not the same haplotype")
                            isUnplaced = False
                            break
                        else:
                            hasUnplaced = True
                    else:
                        print(f"Skipping {row.name} as it is not in scfPathToScaffoldDict")
                        # newList.append([sample,'NA','UNUSED',row.name, row.path])
                        isUnplaced = False
                        break
                    
                
                elif node in mashmapUnplaced:
                    isUnplaced = False
                    if row.name in scfPathToScaffoldDict.keys():
                        isComponent = True
                        # newList.append([sample,scfPathToScaffoldDict[row.name],'UNPLACED_Cmpnt',row.name, row.path])
                        
                    else:
                        print(f"Skipping {row.name} as it is not in scfPathToScaffoldDict")
                        # newList.append([sample,'NA','UNUSED',row.name, row.path])
                        break

                else:
                    if node in utigToChrDict and utigToChrDict[node] != 'chrY':
                        isUnplaced = False
                        if hasUnplaced:
                            print(f"{sample} :: WARNING :: {row.name} has unplaced tigs {listUnplaced} but connected to a non-Y {node} ( {utigToChrDict[node]} ) : {row.path}")
            
                        break

            if hasUnplaced and isUnplaced :
                print(f"Adding {row.name} to unplaced. Path: {row.path}")
                newList.append([sample,scfPathToScaffoldDict[row.name],'UNPLACED',row.name, row.path])
            elif not isUnplaced and isComponent:
                print(f"Adding {row.name} to unplaced component. Path: {row.path}")
                newList.append([sample,scfPathToScaffoldDict[row.name],'UNPLACED_Cmpnt',row.name, row.path])

    
    # Build a dataframe for this newList
    finalDF = pd.DataFrame(data=newList)
    # finalDF.columns=['Genome', 'Seq', 'Selection_Process', 'Path']
    finalDF.to_csv(str(args.output), sep='\t', index=False, header=False)

if __name__=="__main__":
    main()