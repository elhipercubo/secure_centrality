typedef long long ll;
const double infinity = 1e150; 

typedef int vertex;

unordered_map<ll, int> vertices;    // consecutively numbered
vector<ll> vertex_labels;

typedef pair<vertex, vertex> pvv;
typedef vector<pair<pvv, int>> Graph;
//typedef map<pvv, int> Graph;
//typedef vector<vector<vertex>> Graph;  // adjacency lists; adj[i] is sorted in increasing order

template<typename T>
using Matrix = vector<vector<T>>;

/*
 * Returns the integer representing a given vertex label.
 * If there is no integer for this label, assign one to it, in increasing
 * order.
 *
 */
vertex vertex_id(ll label) {
    auto it = vertices.find(label);
    if (it == vertices.end()) {
        int size = vertices.size();
//        adj.push_back(vector<int>());
        vertex_labels.push_back(label);
        return vertices[label] = size;
    } else
        return it->second;
}

/*
 * Returns the next pair of labels in the input in_file
 *
 * Each line in the input file is either a comment if it starts with '#' or a
 * pair of integers, "a b", representing an edge between node 'a' and node 'b'.
 *
 */
bool next_edge(FILE* in_file, array<ll, 3>& e) {
    char line[1000], c;
    // Skip comments
    while ((c = fgetc(in_file)) == '#') {
        fgets(line, sizeof(line), in_file);
    }

    // Handle the case where the last line is a comment and there's no newline.
    if (feof(in_file)) {
        fclose(in_file);
        return false;
    }

    ungetc(c, in_file);

    char *a, *b, *d;
    fgets(line, sizeof(line), in_file);
    a = strtok(line, " \t\r\n");
    b = strtok(NULL, " \t\r\n");
    d = strtok(NULL, " \t\r\n");

    if (a == NULL || b == NULL || strlen(a) == 0 || strlen(b) == 0) {
        if (!feof(in_file)) {
            fprintf(stderr, "Malformed line in input file: %s\n", line);
            exit(1);
        }
        fclose(in_file);
        return false;
    }

    e[0] = atoll(a);
    e[1] = atoll(b);
    if (!d) e[2] = 1;
    else e[2] = atoll(d);
    return true;
}

/*
void read_binary_graph(const char* in_filename, Graph& Graph) {
    FILE* in_file = fopen(in_filename, "r");
    if (!in_file) {
        fprintf(stderr, "Cannot open input file: %s\n", in_filename);
        exit(-1);
    }
    int buf[3], num_edges = 0;
    while (fread(buf, sizeof(int), 3, in_file) == 3) {
        vertex a = vertex_id(buf[0]), b = vertex_id(buf[1]);
        int weight = buf[2];
//        printf("read %i %i %i (%i %i)\n", buf[0], buf[1], weight, a, b);
        Graph[pvv(a, b)] = weight;
        ++num_edges;
    }
    printf("    m = %i\n", num_edges);
}
*/

/*
 * Reads the graph from the input file in_file and fills in the general graph
 * data.
 *
 * Each line in the input file is either a comment if it starts with '#' or a
 * pair of integers, "a b", representing an edge between node 'a' and node 'b'.
 *
 */
void read_graph(const char* in_filename, Graph& graph, bool duplicate_edges) {
    FILE* in_file = fopen(in_filename, "r");
    if (!in_file) {
        fprintf(stderr, "Cannot open input file: %s\n", in_filename);
        exit(-1);
    }

    array<ll, 3> edge;
    int num_edges = 0;
    while (next_edge(in_file, edge)) {
        vertex a = vertex_id(edge[0]), b = vertex_id(edge[1]);
        if (a == b) continue; // we assume no loops, although we could handle them.

        int weight = edge[2];
        graph.push_back(make_pair(pvv(a, b), weight));
        if (duplicate_edges) {
            graph.push_back(make_pair(pvv(b, a), weight));
            ++num_edges;
        }
        ++num_edges;
    }
    if (duplicate_edges) {
        assert(num_edges % 2 == 0);
        num_edges /= 2;
    }
    int n = vertices.size();
    printf("%s: %d edges.\n", in_filename, num_edges);
    printf("    n = %i m = %i avg_deg = %lf\n\n", n, num_edges, 2.0 * num_edges / n);
}
