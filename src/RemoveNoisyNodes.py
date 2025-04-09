import sys
from collections import defaultdict, deque

class Graph:
	def __init__(self):
		self.nodes = {}  # Node ID -> length
		self.edges = defaultdict(list)  # Node+Direction -> list of (target_node+direction, overlap)

	def add_node(self, node_id, length, cov):
		self.nodes[node_id] = [length, cov]
	
	def get_coverage(self, node_id):
		if node_id in self.nodes:
			return self.nodes[node_id][1]
		else:
			print(f":: ERROR :: {node_id} not found in gfa")
			return -1

	def add_edge(self, from_node, from_dir, to_node, to_dir, overlap):
		# Add edges with directionality
		if from_dir == '+':
			from_key = f"{from_node}E"
		else:
			from_key = f"{from_node}B"
		if to_dir == '+':
			to_key = f"{to_node}B"
		else:
			to_key = f"{to_node}E"
		self.edges[from_key].append((to_key, overlap))

	def ignore_node(self, node, avg_cov):
		padding = 13000 # Padding 13 kb: node len - overlap should be > 13kb
		length = self.nodes[node][0] # Get the length of the node
		cov = self.nodes[node][1] # Get the coverage of the node
		cov_threshold = avg_cov / 2 # Set the coverage threshold
		overlaps = {}
		isNoisy = False
		noEdges = 0

		direction = 'B'
		node_key = f"{node}{direction}"
		if node_key in self.edges:
			# print(f"Edges self.edges[{node_key}]: {self.edges[node_key]}")
			for from_node, overlap_from in self.edges[node_key]:
				# keep track of each overlap
				overlaps[from_node] = overlap_from
		else:
			noEdges += 1
		
		direction = 'E'
		node_key = f"{node}{direction}"
		if node_key in self.edges:
			# print(f"Edges self.edges[{node_key}]: {self.edges[node_key]}")
			for to_node, overlap_to in self.edges[node_key]:
				if len(overlaps.items()) > 0:
					# add up the overlaps for each + and - direction
					for from_node, overlap_from in overlaps.items():
						# print(f":: DEBUG :: {node} {cov}x, {length} - ( {from_node} ({overlap_from}) + {to_node} ({overlap_to}) ) = ({length - (overlap_to + overlap_from)}) ")
						print(f"{node} {cov}x, {length} - ( {from_node} ({overlap_from}) + {to_node} ({overlap_to}) ) = ({length - (overlap_to + overlap_from)}) ")
						if cov > 0 and cov < cov_threshold and ( overlap_to + overlap_from + padding > length ):
							isNoisy = True
				else:
					print(f"{node} {cov}x, {length} - {to_node} ({overlap_to}) = ({length - overlap_to}) ")
					if cov > 0 and cov < cov_threshold and ( overlap_to + ( padding / 2 ) > length ):
						isNoisy = True
		else:
			# has only edges from B
			for from_node, overlap_from in overlaps.items():
				print(f"{node} {cov}x, {length} - {from_node} ({overlap_from}) = ({length - overlap_from}) ")
				if cov > 0 and cov < cov_threshold and ( overlap_from + ( padding / 2 ) > length ):
					isNoisy = True
			noEdges += 1
		
		if noEdges == 2:
			# no edges in either direction
			if cov > 0 and cov < cov_threshold and length < padding:
				print(f"{node} {cov}x, {length} <  {padding} - no edges")
				isNoisy = True
			else:
				print(f"{node} {cov}x, {length} >= {padding} - no edges")
		
		return isNoisy

def parse_gfa(file_path):
    graph = Graph()
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if parts[0] == 'S':  # Node line
                node_id = parts[1]
                length = int(parts[3].split(':')[-1])  # Extract length from LN tag
                rc = int(parts[4].split(':')[-1])  # Extract reverse complement from RC tag
                if rc == 0:
                    cov = 0.0
                else:
                    cov = float(rc) / length  # Extract coverage from the tag
                graph.add_node(node_id, length, cov)

            elif parts[0] == 'L':  # Edge line
                from_node = parts[1]
                from_dir = parts[2]
                to_node = parts[3]
                to_dir = parts[4]
                overlap = int(parts[5][:-1])  # Remove trailing 'M' and convert to int
                graph.add_edge(from_node, from_dir, to_node, to_dir, overlap)
    return graph

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("Usage: python RemoveNoisyNodes.py asm.noseq.gfa tigs.csv")
		sys.exit(1)

	gfa_file = sys.argv[1]

	# Parse GFA file
	graph = parse_gfa(gfa_file)
	nodes_file = sys.argv[2]

	# Collect avgerage coverage of all nodes in scaffold path
	with open(nodes_file, 'r') as fi:
		for line in fi:
			tokens = line.strip().split('\t')
			if tokens[2] == 'SCAFFOLD':
				cov_sum = 0
				node_count = 0
				nodes_candidate = tokens[-1].split(',')
				for node in nodes_candidate:
					if node[:4] == 'utig':
						node = node[:-1]
						cov = graph.get_coverage(node)
						if cov > 9 and cov < 100:
							cov_sum += graph.get_coverage(node)
							node_count += 1
	
	if node_count > 0:
		avg_cov = cov_sum / node_count
		print(f"Average coverage of nodes in SCAFFOLD path: {cov_sum} / {node_count} = {avg_cov}")
	else:
		print("No valid nodes found in the SCAFFOLD path")
		sys.exit(1)


	with open(nodes_file, 'r') as fi:
		with open(nodes_file[:-4] + ".no_noise.txt", 'w') as fo:
			for line in fi:
				tokens = line.strip().split('\t')
				if tokens[2] == 'UNPLACED' or tokens[2] == 'UNPLACED_Cmpnt':
					nodes_candidate = tokens[-1].split(',')
					# only include singleton nodes
					if len(nodes_candidate) == 1 and nodes_candidate[0][:4] == 'utig':
						node = nodes_candidate[0][:-1]
						if graph.ignore_node(node, avg_cov):
							print(f"Excluding {line.strip()}")
						else:
							fo.write(line)
					else:
						fo.write(line)
				else:
					fo.write(line)