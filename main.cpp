#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <set>
#include <string>
#include <windows.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <climits>
#include <cctype>

using namespace std;

void clearScreen()
{
    system("CLS");
}

void displayMenu()
{
    cout << "---User Menu---" << endl << endl;
    cout << "1-847. Shortest Path Visiting All Nodes" << endl;
    cout << "2-943. Find the Shortest Superstring" << endl;
    cout << "3-332. Reconstruct Itinerary" << endl;
    cout << "4-2097. Valid Arrangement of Pairs" << endl;
    cout << "5-Harta" << endl;
    cout << "6-Negot" << endl;
    cout << "7-Senat" << endl;
    cout << "8-Paznici" << endl;
    cout << "9-Exit" << endl;
}

class Graph
{
private:

    int n; /// The number of nodes of the graph.
    vector<vector<int>> adjList; /// This vector stores the adjacency list representation of the graph.

    /// Private function to perform BFS with bitmask to find the shortest path
    /// \return The length of the shortest path that visits every node
    int bfsWithBitmask();

    /// Private function to perform depth-first search to reconstruct the itinerary
    /// \param airport The current airport being processed
    /// \param graph The adjacency list representing the flights between airports
    /// \param result The vector to store the reconstructed itinerary
    void dfs(const string& airport, unordered_map<string, multiset<string>>& graph, vector<string>& result);

    /// Private function to perform DFS traversal to find Eulerian path
    /// \param node: Current node being processed
    /// \param parent: Parent node of the current node
    /// \param adjList: Adjacency list representing the graph
    /// \param s: Stack to store the Eulerian path
    void dfsEulerianPath(int node, int parent, unordered_map<int, vector<int>>& adjList, stack<vector<int>>& s);

    /// Private function to perform Breadth-First Search to find an augmenting path from the source to the sink
    /// \param source: Source node
    /// \param sink: Sink node
    /// \param parent: Parent vector to track the augmenting path
    /// \param adjList: Adjacency list representing the graph
    /// \param capacity: Capacities of edges in the graph
    /// \param flow: Flow values in the flow network
    /// \return Returns the minimum flow value in the found augmenting path
    int bfs(int source, int sink, vector<int>& parent, vector<vector<int>>& adjList, vector<vector<int>>& capacity, vector<vector<int>>& flow);

    /// Private function to find the maximum flow in the given flow network using the Ford-Fulkerson algorithm with Edmonds-Karp implementation
    /// \param source: Source node
    /// \param sink: Sink node
    /// \param adjList: Adjacency list representing the graph
    /// \param capacity: Capacities of edges in the graph
    /// \param flow: Flow values in the flow network
    /// \return Returns the maximum flow from source to sink
    int maxFlow(int source, int sink, vector<vector<int>>& adjList, vector<vector<int>>& capacity, vector<vector<int>>& flow);

    /// Private function to perform Breadth-First Search (BFS) to find the minimum cut in a graph
    /// \param source: The source node for the BFS traversal
    /// \param adjList: The adjacency list representation of the graph
    /// \param capacity: The capacity matrix representing the maximum flow that can pass through each edge
    /// \param flow: The flow matrix representing the current flow through each edge
    /// \return Returns a vector of boolean values indicating whether each node is reachable from the source (true if reachable, false otherwise)
    vector<bool> bfsMinCut(int source, vector<vector<int>>& adjList, vector<vector<int>>& capacity, vector<vector<int>>& flow);








public:

    /// Constructor to initialize the graph with given edges and directed flag.
    Graph(int n, vector<vector<int>>& edges, bool directed);
    /// Constructor without parameters.
    Graph();

    /// Public function to find the shortest path that visits every node
    /// \return The length of the shortest path that visits every node
    /// https://leetcode.com/problems/shortest-path-visiting-all-nodes/
    int shortestPathLength();

    /// Public function to return the smallest string that contains each string in words as a substring.
    /// \param words: A vector of strings representing the input words.
    /// \return A string representing the smallest superstring containing all input words.
    /// https://leetcode.com/problems/find-the-shortest-superstring/
    string shortestSuperstring(vector<string>& words);

    /// Public function to reconstruct the itinerary given a list of airline tickets
    /// \param tickets The list of airline tickets represented as pairs of departure and arrival airports
    /// \return The reconstructed itinerary as a vector of strings
    /// https://leetcode.com/problems/reconstruct-itinerary/description/
    vector<string> findItinerary(vector<vector<string>>& tickets);

    /// Public function to find a valid arrangement of pairs using Eulerian path
    /// \param pairs: Input pairs representing the edges of the graph
    /// \return Valid arrangement of pairs as a vector of vectors
    /// https://leetcode.com/problems/valid-arrangement-of-pairs/
    vector<vector<int>> validArrangement(vector<vector<int>>& pairs);

    /// Public function to solve harta problem
    /// \param source: Source node
    /// \param sink: Sink node
    /// \param adjList: Adjacency list representing the graph
    /// \param capacity: Capacities of edges in the graph
    /// \param flow: Flow values in the flow network
    /// \return Returns the maximum flow from source to sink
    /// https://www.infoarena.ro/problema/harta
    int harta(int source, int sink, vector<vector<int>>& adjList, vector<vector<int>>& capacity, vector<vector<int>>& flow);


    /// Public function to solve negot problem
    /// \param source: Source node
    /// \param sink: Sink node
    /// \param adjList: Adjacency list representing the graph
    /// \param capacity: Capacities of edges in the graph
    /// \param flow: Flow values in the flow network
    /// \return Returns the maximum flow from source to sink
    /// https://infoarena.ro/problema/negot
    int negot(int source, int sink, vector<vector<int>>& adjList, vector<vector<int>>& capacity, vector<vector<int>>& flow);

    /// Public function to solve senat problem
    /// \param source: Source node
    /// \param sink: Sink node
    /// \param adjList: Adjacency list representing the graph
    /// \param capacity: Capacities of edges in the graph
    /// \param flow: Flow values in the flow network
    /// \return Returns the maximum flow from source to sink
    /// https://infoarena.ro/problema/senat
    int senat(int source, int sink, vector<vector<int>>& adjList, vector<vector<int>>& capacity, vector<vector<int>>& flow);

    /// Public function to solve paznici problem
    /// \param source: Source node
    /// \param sink: Sink node
    /// \param adjList: Adjacency list representing the graph
    /// \param capacity: Capacities of edges in the graph
    /// \param flow: Flow values in the flow network
    /// \return Returns the maximum flow from source to sink
    /// https://www.infoarena.ro/problema/paznici
    vector<bool> paznici(int source, int sink, vector<vector<int>>& adjList, vector<vector<int>>& capacity, vector<vector<int>>& flow);

};


/// Constructor to initialize the graph with given edges and directed flag.
Graph::Graph(int n, vector<vector<int>>& edges, bool directed)
{
    this->n = n; /// Initialize the number of nodes.
    this->adjList.resize(this->n + 1); /// Resize the adjacency list to accommodate 'n' nodes.

    /// Iterate through the given edges and populate the adjacency list.
    for (const auto& edge : edges)
    {
        this->adjList[edge[0]].push_back(edge[1]); /// Add edge from node 'edge[0]' to node 'edge[1]'.

        /// If the graph is undirected, add reverse edge from 'edge[1]' to 'edge[0]'.
        if (directed == false)
        {
            this->adjList[edge[1]].push_back(edge[0]);
        }
    }
}

/// Constructor without parameters.
Graph::Graph()
{
    this->n = 0; /// Initialize the number of nodes to 0 for the default constructor.
    this->adjList = vector<vector<int>>(); /// Initialize the adjacency list as an empty vector of vectors.
}

/// Private function to perform BFS with bitmask to find the shortest path
/// \return The length of the shortest path that visits every node
int Graph::bfsWithBitmask()
{
    /// Calculate the bitmask for the target state where all nodes are visited
    int targetMask = (1 << this->n) - 1;

    /// Queue for BFS
    queue<pair<int, int>> q;

    /// 2D vector to track visited nodes with specific bitmask
    vector<vector<int>> visited(this->n, vector<int>(1 << this->n, false));

    /// Initialize the queue and mark the starting nodes as visited
    for (int i = 0; i < this->n; i++)
    {
        q.push({i, 1 << i});
        visited[i][1 << i] = true;
    }

    /// Counter for the number of steps
    int cnt = 0;

    /// Perform BFS
    while (!q.empty())
    {
        /// Process nodes at the current level
        int size = q.size();
        for (int i = 0; i < size; i++)
        {
            /// Get the current node and bitmask
            int node = q.front().first;
            int mask = q.front().second;
            q.pop();

            /// Check if the target state is reached
            if (mask == targetMask)
            {
                return cnt;
            }

            /// Explore neighbors and update the bitmask
            for (int neighbor : this->adjList[node])
            {
                int nextMask = mask | (1 << neighbor);

                /// If the neighbor with the updated bitmask is not visited, add to the queue
                if (!visited[neighbor][nextMask])
                {
                    q.push({neighbor, nextMask});
                    visited[neighbor][nextMask] = true;
                }
            }
        }

        /// Move to the next level
        cnt++;
    }

    /// If the target state is not reached, return -1
    return -1;
}


/// Private function to perform depth-first search to reconstruct the itinerary
/// \param airport The current airport being processed
/// \param graph The adjacency list representing the flights between airports
/// \param result The vector to store the reconstructed itinerary
void Graph::dfs(const string& airport, unordered_map<string, multiset<string>>& graph, vector<string>& result)
{
    /// Continue the DFS while there are available flights from the current airport
    while (!graph[airport].empty())
    {
        /// Get the next airport in lexical order
        string nextAirport = *(graph[airport].begin());
        /// Remove the chosen flight from the multiset
        graph[airport].erase(graph[airport].begin());
        /// Recursively perform DFS on the next airport
        this->dfs(nextAirport, graph, result);
    }

    /// Add the current airport to the result vector
    result.push_back(airport);
}

/// Private function to perform DFS traversal to find Eulerian path
/// \param node: Current node being processed
/// \param parent: Parent node of the current node
/// \param adjList: Adjacency list representing the graph
/// \param s: Stack to store the Eulerian path
void Graph::dfsEulerianPath(int node, int parent, unordered_map<int, vector<int>>& adjList, stack<vector<int>>& s)
{
    /// Traverse outgoing edges from the current node
    while (!adjList[node].empty())
    {
        /// Get the next node
        int nextNode = adjList[node].back();
        /// Remove the chosen edge
        adjList[node].pop_back();
        /// Recursively perform DFS on the next node
        dfsEulerianPath(nextNode, node, adjList, s);
    }

    /// Add the current edge to the Eulerian path
    if (parent != -1)
        s.push({parent, node});
}

/// Private function to perform Breadth-First Search to find an augmenting path from the source to the sink
/// \param source: Source node
/// \param sink: Sink node
/// \param parent: Parent vector to track the augmenting path
/// \param adjList: Adjacency list representing the graph
/// \param capacity: Capacities of edges in the graph
/// \param flow: Flow values in the flow network
/// \return Returns the minimum flow value in the found augmenting path
int Graph::bfs(int source, int sink, vector<int>& parent, vector<vector<int>>& adjList, vector<vector<int>>& capacity, vector<vector<int>>& flow)
{
    /// Initialize parent vector to track the augmenting path
    fill(parent.begin(), parent.end(), -1);
    parent[source] = source;

    /// Use a queue for BFS traversal with pair (node, flow)
    queue<pair<int, int>> q;
    q.push({source, INT_MAX});

    while (!q.empty())
    {
        int currentNode = q.front().first;
        int currentFlow = q.front().second;
        q.pop();

        /// Explore neighbors of the current node
        for (int neighbor : adjList[currentNode])
        {
            /// Check if the neighbor is unvisited and there is available capacity
            if (parent[neighbor] == -1 && capacity[currentNode][neighbor] > flow[currentNode][neighbor])
            {
                /// Update parent and find the minimum flow in the path
                parent[neighbor] = currentNode;
                int newFlow = min(currentFlow, capacity[currentNode][neighbor] - flow[currentNode][neighbor]);

                /// If the sink is reached, return the minimum flow value
                if (neighbor == sink)
                {
                    return newFlow;
                }

                /// Enqueue the neighbor with the new flow value
                q.push({neighbor, newFlow});
            }
        }
    }

    /// If no augmenting path is found, return 0
    return 0;
}

/// Private function to find the maximum flow in the given flow network using the Ford-Fulkerson algorithm with Edmonds-Karp implementation
/// \param source: Source node
/// \param sink: Sink node
/// \param adjList: Adjacency list representing the graph
/// \param capacity: Capacities of edges in the graph
/// \param flow: Flow values in the flow network
/// \return Returns the maximum flow from source to sink
int Graph::maxFlow(int source, int sink, vector<vector<int>>& adjList, vector<vector<int>>& capacity, vector<vector<int>>& flow)
{
    /// Initialize parent vector for augmenting paths
    vector<int> parent(capacity.size());

    /// Initialize maximum flow to 0
    int maxFlow = 0;
    int newFlow;

    /// Repeat while there is an augmenting path
    while (newFlow = this->bfs(source, sink, parent, adjList, capacity, flow))
    {
        /// Update the maximum flow with the new flow value
        maxFlow = maxFlow + newFlow;

        int currentNode = sink;
        /// Update flow values in the augmenting path
        while (currentNode != source)
        {
            int previousNode = parent[currentNode];
            flow[previousNode][currentNode] = flow[previousNode][currentNode] + newFlow;
            flow[currentNode][previousNode] = flow[currentNode][previousNode] - newFlow;
            currentNode = previousNode;
        }
    }

    /// Return the final maximum flow
    return maxFlow;
}

/// Private function to perform Breadth-First Search (BFS) to find the minimum cut in a graph
/// \param source: The source node for the BFS traversal
/// \param adjList: The adjacency list representation of the graph
/// \param capacity: The capacity matrix representing the maximum flow that can pass through each edge
/// \param flow: The flow matrix representing the current flow through each edge
/// \return Returns a vector of boolean values indicating whether each node is reachable from the source (true if reachable, false otherwise)
vector<bool> Graph::bfsMinCut(int source, vector<vector<int>>& adjList, vector<vector<int>>& capacity, vector<vector<int>>& flow)
{
    /// Initialize a vector to keep track of visited nodes during BFS
    vector<bool> visited(capacity.size(), false);

    /// Mark the source node as visited and enqueue it
    visited[source] = true;
    queue<int> q;
    q.push(source);

    /// Perform BFS
    while (!q.empty())
    {
        /// Dequeue a node and process its neighbors
        int currentNode = q.front();
        q.pop();

        /// Iterate over the neighbors of the current node
        for (int neighbor : adjList[currentNode])
        {
            /// Check if the neighbor is not visited
            if (visited[neighbor] == false && capacity[currentNode][neighbor] > flow[currentNode][neighbor])
            {
                /// Mark the neighbor as visited and enqueue it for further exploration
                visited[neighbor] = true;
                q.push(neighbor);
            }
        }
    }

    /// Return the vector indicating reachability from the source
    return visited;
}



















/// Public function to find the shortest path that visits every node
/// \return The length of the shortest path that visits every node
/// https://leetcode.com/problems/shortest-path-visiting-all-nodes/
int Graph::shortestPathLength()
{
    /// If there is only one node, the shortest path is 0
    if (this->n == 1)
        return 0;

    /// Call the BFS function to find the shortest path
    return this->bfsWithBitmask();
}


/// Public function to return the smallest string that contains each string in words as a substring.
/// \param words: A vector of strings representing the input words.
/// \return A string representing the smallest superstring containing all input words.
/// https://leetcode.com/problems/find-the-shortest-superstring/
string Graph::shortestSuperstring(vector<string>& words)
{
    /// Get the number of words in the input
    int n = words.size();

    /// Matrix to store the overlap lengths between words
    vector<vector<int>> overlap(n, vector<int>(n, 0));

    /// Calculate overlap lengths for each pair of words
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                /// Find the minimum length for comparison
                int len = min(words[i].size(), words[j].size());

                /// Check for overlap by comparing suffix and prefix
                for (int k = len; k >= 0; k--)
                {
                    if (words[i].substr(words[i].size() - k) == words[j].substr(0, k))
                    {
                        overlap[i][j] = k;  /// Store the overlap length
                        break;
                    }
                }
            }
        }
    }

    /// Dynamic programming matrix to store minimum superstring lengths
    vector<vector<int>> dp(1 << n, vector<int>(n, INT_MAX));

    /// Matrix to store parent indices for reconstruction
    vector<vector<int>> parent(1 << n, vector<int>(n, -1));

    /// Initialize base cases for dynamic programming
    for (int i = 0; i < n; ++i)
    {
        dp[1 << i][i] = words[i].size();
    }

    /// Dynamic programming to find minimum superstring lengths
    for (int mask = 1; mask < (1 << n); mask++)
    {
        for (int i = 0; i < n; i++)
        {
            /// Check if word i is in the current superstring represented by the mask
            if ((mask & (1 << i)) != 0)
            {
                for (int j = 0; j < n; j++)
                {
                    /// Check if word j is not yet in the current superstring
                    if ((mask & (1 << j)) == 0)
                    {
                        /// Calculate new length for adding word j to the current superstring
                        int newLen = dp[mask][i] + words[j].size() - overlap[i][j];

                        /// Update DP matrix if the new length is smaller
                        if (newLen < dp[mask | (1 << j)][j])
                        {
                            dp[mask | (1 << j)][j] = newLen;
                            parent[mask | (1 << j)][j] = i;
                        }
                    }
                }
            }
        }
    }

    /// Find the minimum length and corresponding ending node for reconstruction
    int minLen = INT_MAX;
    int endNode = -1;

    for (int i = 0; i < n; i++)
    {
        if (dp[(1 << n) - 1][i] < minLen)
        {
            minLen = dp[(1 << n) - 1][i];
            endNode = i;
        }
    }

    /// Reconstruct the superstring using parent indices
    string result = words[endNode];

    int mask = (1 << n) - 1;
    while (parent[mask][endNode] != -1)
    {
        int prevNode = parent[mask][endNode];
        /// Append the necessary portion of the previous word to the result
        result = words[prevNode].substr(0, words[prevNode].size() - overlap[prevNode][endNode]) + result;
        /// Update mask and ending node for the next iteration
        mask = mask ^ (1 << endNode);
        endNode = prevNode;
    }

    /// Return the final reconstructed superstring
    return result;
}


/// Public function to reconstruct the itinerary given a list of airline tickets
/// \param tickets The list of airline tickets represented as pairs of departure and arrival airports
/// \return The reconstructed itinerary as a vector of strings
/// https://leetcode.com/problems/reconstruct-itinerary/description/
vector<string> Graph::findItinerary(vector<vector<string>>& tickets)
{
    /// Adjacency list to represent the flights between airports
    unordered_map<string, multiset<string>> graph;
    /// Vector to store the reconstructed itinerary
    vector<string> result;

    /// Populate the adjacency list based on the provided tickets
    for (const auto& ticket : tickets)
    {
        graph[ticket[0]].insert(ticket[1]);
    }

    /// Start the DFS from the "JFK" airport
    this->dfs("JFK", graph, result);

    /// Reverse the result vector to obtain the correct order
    reverse(result.begin(), result.end());

    /// Return the reconstructed itinerary
    return result;
}

/// Public function to find a valid arrangement of pairs using Eulerian path
/// \param pairs: Input pairs representing the edges of the graph
/// \return Valid arrangement of pairs as a vector of vectors
/// https://leetcode.com/problems/valid-arrangement-of-pairs/
vector<vector<int>> Graph::validArrangement(vector<vector<int>>& pairs)
{
    /// Adjacency list for the graph
    unordered_map<int, vector<int>> adjList;

    /// Degree of each node
    unordered_map<int, int> degree;

    /// Starting node for Eulerian path
    int startingNode = -1;

    /// Populate adjacency list and degree based on input pairs
    for (auto pair : pairs)
    {
        adjList[pair[0]].push_back(pair[1]);
        degree[pair[0]]++;
        degree[pair[1]]--;
    }

    /// Find a node with positive outdegree as the starting node
    for (auto nodeDegree : degree)
    {
        if (nodeDegree.second > 0)
        {
            startingNode = nodeDegree.first;
            break;
        }
    }

    /// If no such node found, arbitrarily choose the starting node
    if (startingNode == -1)
        startingNode = pairs[0][0];

    /// Stack for DFS traversal
    stack<vector<int>> s;

    /// Perform DFS to find Eulerian path
    dfsEulerianPath(startingNode, -1, adjList, s);

    /// Result vector to store the valid arrangement
    vector<vector<int>> result;

    /// Transfer elements from stack to result vector
    while (!s.empty())
    {
        result.push_back(s.top());
        s.pop();
    }

    return result;
}

/// Public function to solve harta problem
/// \param source: Source node
/// \param sink: Sink node
/// \param adjList: Adjacency list representing the graph
/// \param capacity: Capacities of edges in the graph
/// \param flow: Flow values in the flow network
/// \return Returns the maximum flow from source to sink
/// https://www.infoarena.ro/problema/harta
int Graph::harta(int source, int sink, vector<vector<int>>& adjList, vector<vector<int>>& capacity, vector<vector<int>>& flow)
{
    return this->maxFlow(source, sink, adjList, capacity, flow);
}

/// Public function to solve negot problem
/// \param source: Source node
/// \param sink: Sink node
/// \param adjList: Adjacency list representing the graph
/// \param capacity: Capacities of edges in the graph
/// \param flow: Flow values in the flow network
/// \return Returns the maximum flow from source to sink
/// https://infoarena.ro/problema/negot
int Graph::negot(int source, int sink, vector<vector<int>>& adjList, vector<vector<int>>& capacity, vector<vector<int>>& flow)
{
    return this->maxFlow(source, sink, adjList, capacity, flow);
}

/// Public function to solve senat problem
/// \param source: Source node
/// \param sink: Sink node
/// \param adjList: Adjacency list representing the graph
/// \param capacity: Capacities of edges in the graph
/// \param flow: Flow values in the flow network
/// \return Returns the maximum flow from source to sink
/// https://www.infoarena.ro/problema/senat
int Graph::senat(int source, int sink, vector<vector<int>>& adjList, vector<vector<int>>& capacity, vector<vector<int>>& flow)
{
    return this->maxFlow(source, sink, adjList, capacity, flow);
}

/// Public function to solve paznici problem
/// \param source: Source node
/// \param sink: Sink node
/// \param adjList: Adjacency list representing the graph
/// \param capacity: Capacities of edges in the graph
/// \param flow: Flow values in the flow network
/// \return Returns the maximum flow from source to sink
/// https://www.infoarena.ro/problema/paznici
vector<bool> Graph::paznici(int source, int sink, vector<vector<int>>& adjList, vector<vector<int>>& capacity, vector<vector<int>>& flow)
{
    this->maxFlow(source, sink, adjList, capacity, flow);
    return this->bfsMinCut(source,adjList,capacity,flow);
}




int main()
{
    int cnt=0;
    while(true)
    {
        displayMenu();
        int command;
        cin>>command;
        switch(command)
        {
        case 1:
        {
            clearScreen();

            cout<<"n= ";
            int n;
            cin>>n;

            int x, y;
            cout<<endl<<"edges (enter -1 to stop)= ";
            vector<vector<int>> edges;
            while (cin>>x && x!=-1 && cin>>y)
            {
                edges.push_back({x, y});
            }
            cout<<endl;


            Graph G(n,edges,false);

            cout<<G.shortestPathLength();

            cout<<endl;
            cout<<"Press 'Enter' to return to the menu."<<endl;
            cin.ignore();
            while(cin.get() != '\n');
            clearScreen();
            break;


        }

        case 2:
        {
            clearScreen();

            cout << "words (enter -1 to stop)= ";
            string word;

            vector<string> words;

            while (cin >> word && word != "-1")
            {
                words.push_back(word);
            }
            cout<<endl;

            Graph G;
            cout<<G.shortestSuperstring(words);


            cout<<endl;
            cout<<"Press 'Enter' to return to the menu."<<endl;
            cin.ignore();
            while(cin.get() != '\n');
            clearScreen();
            break;
        }

        case 3:
        {
            clearScreen();

            cout << "tickets (enter -1 to stop)= ";
            string x,y;

            vector<vector<string>> tickets;

            while (cin >> x && x != "-1" && cin >> y)
            {
                tickets.push_back({x,y});
            }
            cout<<endl;

            Graph G;
            vector<string> result;
            result=G.findItinerary(tickets);

            for(int i=0; i < result.size(); i++)
                cout<<result[i] << " ";

            cout<<endl;
            cout<<"Press 'Enter' to return to the menu."<<endl;
            cin.ignore();
            while(cin.get() != '\n');

            clearScreen();
            break;
        }
        case 4:
        {
            clearScreen();

            int x, y;
            cout<<"pairs (enter -1 to stop)= ";
            vector<vector<int>> pairs;
            while (cin>>x && x!=-1 && cin>>y)
            {
                pairs.push_back({x, y});
            }
            cout<<endl;


            Graph G;
            vector<vector<int>> result;
            result=G.validArrangement(pairs);

            for(int i=0; i < result.size(); i++)
            {
                for(int j=0; j < result[i].size(); j++)
                    cout<<result[i][j] << " ";
                cout<<endl;
            }


            cout<<endl;
            cout<<"Press 'Enter' to return to the menu."<<endl;
            cin.ignore();
            while(cin.get() != '\n');
            clearScreen();
            break;
        }
        case 5:
        {
            clearScreen();


            int n;
            cout<<"n= ";
            cin >> n;
            cout<<endl;

            int source = 0;
            int sink = n * 2 + 1;

            vector<vector<int>> adjList(n * 2 + 2, vector<int>(n * 2 + 2, 0));
            vector<vector<int>> capacity(n * 2 + 2, vector<int>(n * 2 + 2, 0));
            vector<vector<int>> flow(n * 2 + 2, vector<int>(n * 2 + 2, 0));

            cout<<"x, y= ";

            for (int i = 1; i <= n; i++)
            {
                int outDegree, inDegree;
                cin >> outDegree >> inDegree;

                adjList[source].push_back(i);
                adjList[i].push_back(source);
                capacity[source][i] = outDegree;

                adjList[sink].push_back(i + n);
                adjList[i+n].push_back(sink);
                capacity[i + n][sink] = inDegree;

            }

            for(int i = 1; i <= n ; i++)
                for(int j = 1 ; j <= n ; j++)
                    if(i != j)
                    {
                        adjList[i].push_back(j + n);
                        adjList[j+n].push_back(i);
                        capacity[i][j+n] = 1;
                    }


            cout<<endl;
            Graph G;
            cout<<G.harta(source,sink,adjList,capacity,flow)<<endl;

            for(int i = 1 ; i <= n ; i++)
                for(int j = n + 1 ; j <= n+n ; j++)
                    if(flow[i][j])
                        cout<<i<<" "<<j-n<<endl;





            cout<<"Press 'Enter' to return to the menu."<<endl;
            cin.ignore();
            while(cin.get() != '\n');
            clearScreen();
            break;
        }
        case 6:
        {
            clearScreen();



            cout<<"n, m, k= ";
            int n, m, k;
            cin >> n >> m >> k;

            int source = 0;
            int sink = n + m + 1;

            vector<vector<int>> adjList(n + m + 2);
            vector<vector<int>> capacity(n + m + 2, vector<int>(n + m + 2, 0));
            vector<vector<int>> flow(n + m + 2, vector<int>(n + m + 2, 0));

            cout<<endl;


            for (int i = 1; i <= n; i++)
            {
                cout<<"k= ";
                int x;
                cin >> x;

                adjList[source].emplace_back(i);
                adjList[i].emplace_back(source);
                capacity[source][i] = k;

                for (int j = 0; j < x; j++)
                {
                    cout<<"T"<<j<<"= ";
                    int y;
                    cin >> y;
                    adjList[i].emplace_back(y + n);
                    adjList[y + n].emplace_back(i);
                    capacity[i][y + n] = 1;
                }
            }

            for (int i = n + 1; i <= n + m; i++)
            {
                adjList[sink].emplace_back(i);
                adjList[i].emplace_back(sink);
                capacity[i][sink] = 1;
            }

            cout<<endl;
            Graph G;
            cout << G.negot(source, sink, adjList, capacity, flow);





            cout<<endl;
            cout<<"Press 'Enter' to return to the menu."<<endl;
            cin.ignore();
            while(cin.get() != '\n');
            clearScreen();
            break;
        }
        case 7:
        {
            clearScreen();

            cout<<"n, m= ";
            int n, m;
            cin >> n >> m;
            cout<<endl;

            int source = 0;
            int sink = n + m + 1;

            vector<vector<int>> adjList(n + m + 2);
            vector<vector<int>> capacity(n + m + 2, vector<int>(n + m + 2, 0));
            vector<vector<int>> flow(n + m + 2, vector<int>(n + m + 2, 0));

            for (int i = 1; i <= n; i++)
            {
                adjList[source].push_back(i);
                adjList[i].push_back(source);
                capacity[source][i] = 1;
            }

            for (int i = n + 1; i <= n + m; i++)
            {
                adjList[i].push_back(sink);
                adjList[sink].push_back(i);
                capacity[i][sink] = 1;
            }

            cout<<"m lines= ";
            string line;
            getline(cin,line);

            for (int i = 1; i <= m; i++)
            {
                getline(cin,line);
                int len = line.size();
                int x = 0;
                for (int j = 0; j < len; j++)
                    if (isdigit(line[j]))
                    {
                        x = x*10 + (line[j] - '0');
                    }
                    else
                    {
                        adjList[x].push_back(n+i);
                        adjList[n+i].push_back(x);
                        capacity[x][n+i] = 1;
                        x = 0;
                    }
                if (x!=0)
                {
                    adjList[x].push_back(n+i);
                    adjList[n+i].push_back(x);
                    capacity[x][n+i] = 1;
                }
            }
            cout<<endl;

            Graph G;
            int finalFlow = G.senat(source,sink,adjList,capacity,flow);

            if (finalFlow != m)
                cout<<0;
            else
            {
                for(int i = n+1 ; i <= n+m ; i++)
                {
                    for(int j = 1 ; j <= n + 1 ; j++)
                        if(flow[j][i])
                            cout<<j<<endl;
                }
            }



            cout<<"Press 'Enter' to return to the menu."<<endl;
            cin.ignore();
            while(cin.get() != '\n');

            clearScreen();
            break;
        }
        case 8:
        {
            clearScreen();



            cout<<"n, m= ";
            int n, m;
            cin >> n >> m;

            int source = 0;
            int sink = n + m + 1;

            vector<vector<int>> adjList(n + m + 2);
            vector<vector<int>> capacity(n + m + 2, vector<int>(n + m + 2, 0));
            vector<vector<int>> flow(n + m + 2, vector<int>(n + m + 2, 0));

            cout<<endl;

            cout<<"n lines= ";

            for (int i = 1; i <= n; i++)
            {
                adjList[source].push_back(i);
                adjList[i].push_back(source);
                capacity[source][i] = 1;
            }

            for (int i = n + 1; i <= n + m; i++)
            {
                adjList[i].push_back(sink);
                adjList[sink].push_back(i);
                capacity[i][sink] = 1;
            }

            string line;
            getline(cin,line);

            for(int i=1; i<=n; i++)
            {
                getline(cin,line);
                for (int j=1; j<=m; j++)
                {
                    if(line[j-1]=='1')
                    {
                        adjList[i].push_back(n+j);
                        adjList[n+j].push_back(i);
                        capacity[i][n+j]=1;
                    }
                }
            }

            cout<<endl;
            Graph G;

            vector<bool> result = G.paznici(source,sink,adjList,capacity,flow);

            for(int i=0; i<n+m+1; i++)
            {
                if(result[i+1]==false && i<n)
                    cout<<char('A'+i);
                if(result[i]==true && i>=n+1)
                    cout<<char('a'+i-n-1);


            }




            cout<<endl;
            cout<<"Press 'Enter' to return to the menu."<<endl;
            cin.ignore();
            while(cin.get() != '\n');

            clearScreen();
            break;
        }

        case 9:
        {
            clearScreen();
            cnt=1;
            break;

        }

        }
        if(cnt==1)
        {
            clearScreen();
            break;
        }
    }

    return 0;
}
