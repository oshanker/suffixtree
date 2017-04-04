package suffix;

import java.io.*;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Iterator;

/**
 * This is a java version based on Mark Nelson's C++ implementation of Ukkonnnen's
 * algorithm (http://marknelson.us/1996/08/01/suffix-trees/).
 * @author oshanker
 *
 */
class SuffixTree {

	private byte[] T ;
	
	
	/**
	 * The input buffer and character count. Please note that N is the length of
	 * the input string -1, which means it denotes the maximum index in the
	 * input buffer.
	 **/
	private int N;
	
	/**
	 * The array of defined nodes. 
	 **/
	private Node[] nodes;
	private Suffix active;
	private static byte endIndex;
	static char endChar;
	static char[] alphabet;
	public static int Count = 0;
	int[] leafCount;
	int treeHeight = 0;
	NumberFormat intFormat = NumberFormat.getIntegerInstance();

	/**
	 * The  node contains the suffix link. Each suffix in
	 * the tree that ends at a particular node can find the next smaller suffix by
	 * following the suffix_node link to a new node. Nodes are stored in a simple
	 * array.
	 **/
	class Node implements Iterable<Edge>{
		class EdgeIterator implements Iterator<Edge>{
			int currentIndex = 0;
			int maxIndex;

			public EdgeIterator() {
				this.currentIndex = 0;
				for(maxIndex=edges.length-1;maxIndex>0;maxIndex--){
					if(edges[maxIndex] != null){
						break;
					}
				}
				if(edges[maxIndex] == null){
					maxIndex = -1;
				}
			}

			public boolean hasNext() {
				return currentIndex <= maxIndex;
			}

			public Edge next() {
				for (; currentIndex  <= maxIndex; currentIndex++) {
					if(edges[currentIndex] != null){
						return edges[currentIndex++];
					}
				}
				return null;
			}
		}
		
		public int suffix_node;
		public int implicitLeafCount;
		Edge[] edges = new Edge[alphabet.length];

		/** Constructor **/
		public Node() {
			suffix_node = -1;
		}
		public Edge find( byte c) {
			if(c==alphabet.length){return null;}
			return edges[c];
		}
		void insert(Edge edge){
			byte idx = T[edge.first_char_index];
			edges[idx] = edge;
		}
		public Iterator<Edge> iterator() {
			return new EdgeIterator();
		}
	}

	/**
	 * When a new tree is added to the table, we step through all the currently
	 * defined suffixes from the active point to the end point. This structure
	 * defines a Suffix by its final character. In the canonical representation,
	 * we define that last character by starting at a node in the tree, and
	 * following a string of characters, represented by first_char_index and
	 * last_char_index. The two indices point into the input string. Note that
	 * if a suffix ends at a node, there are no additional characters needed to
	 * characterize its last character position. When this is the case, we say
	 * the node is Explicit, and set first_char_index > last_char_index to flag
	 * that.
	 */
	class Suffix {
		public int origin_node, first_char_index, last_char_index;

		public Suffix(int node, int start, int stop) {
			origin_node = node;
			first_char_index = start;
			last_char_index = stop;
		}

		public boolean explicit() {
			return first_char_index > last_char_index;
		}

		/**
		 * A suffix in the tree is denoted by a Suffix structure
		 * that denotes its last character. The canonical
		 * representation of a suffix for this algorithm requires
		 * that the origin_node by the closest node to the end
		 * of the tree. To force this to be true, we have to
		 * slide down every edge in our current path until we
		 * reach the final node
		 **/
		public void canonize() {
			if (!explicit()) {
				Edge edge = nodes[origin_node].find( T[first_char_index]);
				int edge_span = edge.last_char_index - edge.first_char_index;
				while (edge_span <= (last_char_index - first_char_index)) {
					first_char_index = first_char_index + edge_span + 1;
					origin_node = edge.end_node;
					if (first_char_index <= last_char_index) {
						edge = nodes[edge.end_node].find( T[first_char_index]);
						edge_span = edge.last_char_index - edge.first_char_index;
					}
				}
			}
		}
	}

	/**
	 * The suffix tree is made up of edges connecting nodes. Each edge
	 * represents a string of characters starting at first_char_index and ending
	 * at last_char_index. 
	 * 
	 * Class Edge
	 **/
	class Edge {
		@Override
		public String toString() {
			return "Edge [first_char_index=" + first_char_index + " " + T[first_char_index] +", last_char_index=" + last_char_index + ", end_node="
					+ end_node + ", start_node=" + start_node + "]";
		}


		public int first_char_index;
		public int last_char_index;
		public int end_node;
		public int start_node;

		/**
		 * I create new edges in the program while walking up the set of
		 * suffixes from the active point to the endpoint. Each time I create a
		 * new edge, I also add a new node for its end point. 
		 **/
		public Edge() {
			start_node = -1;
		}

		/** Constructor **/
		public Edge(int init_first, int init_last, int parent_node) {
			this(init_first, init_last);
			start_node = parent_node;
			if (last_char_index<N){
			   nodes[Count] = new Node();
			}
			end_node = Count++;
			nodes[parent_node].insert(this);
		}

		public Edge(int init_first, int init_last) {
			first_char_index = init_first;
			last_char_index = init_last;
			//start_node = parent_node;
		}
		
		Edge newEdgeWithSpecifiedEndNode(int init_first, int init_last, int child_node){
			Edge newEdge = new Edge(init_first, init_last);
			nodes[Count] = new Node();
			newEdge.start_node = Count++;
			newEdge.end_node = child_node;
			nodes[newEdge.start_node].insert(newEdge);
			return newEdge;
		}


		/**
		 * function SplitEdge ()
		 * When a suffix ends on an implicit node, adding a new character means
		 * I have to split an existing edge. This function is called to split an
		 * edge at the point defined by the Suffix argument. The existing edge
		 * loses its child, as well as some of its trailing characters. The
		 * newly created edge descends to the original child, and now has the
		 * existing edge as a parent. 
		 * The number of characters kept in the original edge and given to
		 * the new node is equal to the number of characters in the suffix
		 * argument, which is last - first + 1. The rest is given to the new edge.
		 **/
		public int splitEdge(Suffix s) {
			if(s.origin_node != start_node){ 
				throw new IllegalStateException(" s.origin_node !=start_node " + s.origin_node + "!=" + start_node );
			}
			//from newnode to end,lastchar 
			Edge second = newEdgeWithSpecifiedEndNode(first_char_index + s.last_char_index - s.first_char_index + 1, 
					last_char_index, end_node);
			//from st, firstchar to newnode
			last_char_index = first_char_index + s.last_char_index - s.first_char_index;
			int newnode = second.start_node;
			end_node = newnode;
			nodes[newnode].suffix_node = s.origin_node;
			
			return newnode;
		}

	}

	/** Constructor 
	 * @param str */
	public SuffixTree(byte[] str) {
		intFormat.setMinimumIntegerDigits(7);
		intFormat.setGroupingUsed(false);
		nodes = new Node[str.length * 2];
		T = str;
		/**
		 * The active point is the first non-leaf
		 * suffix in the tree. We start by setting
		 * this to be the empty string at node 0.
		 * The AddPrefix() function will update this
		 * value after every new prefix is added.
		 **/
		active = new Suffix(0, 0, -1);
		this.N = this.T.length - 1;
		nodes[Count++] = new Node();
	}

	/**
	 * This fills in the leaf counts for each node.
	 * @param start_node
	 * @param level
	 * @return
	 */
	int makeTreeExplicit(int start_node,  int level) {
		int edges = 0;
		int leafCountForNode = nodes[start_node].implicitLeafCount;
		Node current = nodes[start_node];
		for (Edge edge : current) {
				edges++;
				if (edge.last_char_index>=N){
					// this edge points to a leaf. No need to make the recursive call,
					// just handle the logic from within this call
					leafCountForNode++;
					leafCount[edge.end_node] = 1;
					if(level+1>treeHeight){treeHeight= level+1;}
					continue;
				} 
				int childLeafCount = makeTreeExplicit(edge.end_node, level+1);
				leafCountForNode += childLeafCount;
				if (childLeafCount == 1) {
					//returned from a leaf node
					throw new IllegalStateException(" Unexpected leaf node " + edge.end_node);
				}
		}
		if (edges == 0) {
			throw new IllegalStateException("unexpected leaf " + start_node );
		} 
		leafCount[start_node] = leafCountForNode;
		if(level>treeHeight){treeHeight= level;}

		return leafCountForNode;
	}

	/**
	 * This routine constitutes the heart of the algorithm. It is called
	 * repetitively, once for each of the prefixes of the input string. The
	 * prefix in question is denoted by the index of its last character. At each
	 * prefix, we start at the active point, and add a new edge denoting the new
	 * last character, until we reach a point where the new edge is not needed
	 * due to the presence of an existing edge starting with the new last
	 * character. This point is the end point. Luckily for us, the end point
	 * just happens to be the active point for the next pass through the tree.
	 * All we have to do is update it's last_char_index to indicate that it has
	 * grown by a single character, and then this routine can do all its work
	 * one more time. Note that I never create nodes explicitly for leaf nodes.
	 * This is because in a suffix tree, once a leaf, always a leaf. It has no
	 * information which needs to be sortored.
	 **/
	public  void addPrefix(Suffix active, int last_char_index) {
		int parent_node;
		int last_parent_node = -1;
		byte c = last_char_index<T.length?T[last_char_index]:(byte)alphabet.length;
		if(c == alphabet.length){
			if(last_char_index < N){
				throw new IllegalArgumentException(last_char_index + " is not last idx " );
			}
		}
		for (;;) {
			Edge edge;
			parent_node = active.origin_node; 
			if(c==alphabet.length){
				if(parent_node==0 && nodes[parent_node].implicitLeafCount>0){
					break;
				} 
			}
			
			/**
			 * Step 1 is to try and find a matching edge for the given
			 * node. If a matching edge
			 * exists, we are done adding
			 * edges, so we break out of
			 * this big loop.
			 **/
			if (active.explicit()) {
				edge = nodes[active.origin_node].find( c);
				if (edge != null){
					break;
				}
			} else {
				edge = nodes[active.origin_node].find( T[active.first_char_index]);
				int span = active.last_char_index - active.first_char_index;
				if (T[edge.first_char_index + span + 1] == c){
					// a no-op
					break;
				}
				parent_node = edge.splitEdge(active);
			} 
			/**
			 * We didn't find a matching edge, so we create a new one, add
			 * it to the tree at the parent node position, and insert it
			 * into the hash table. When we create a new node, it also means
			 * we need to create a suffix link to the new node from the last
			 * node we visited.
			 **/
			if(c==alphabet.length){
				nodes[parent_node].implicitLeafCount++;
			} else {
				Edge new_edge = new Edge(last_char_index, N, parent_node);
				if (last_parent_node > 0) {
					nodes[last_parent_node].suffix_node = parent_node;
				}
				last_parent_node = parent_node;
			}
			/**
			 * This final step is where we move to the next smaller suffix
			 * using suffix link
			 **/
			if (active.origin_node == 0) {
				active.first_char_index++;
			} else {
				active.origin_node = nodes[active.origin_node].suffix_node;
			}
			active.canonize();
		}
		if(c!=alphabet.length){
			if (last_parent_node > 0) {
				nodes[last_parent_node].suffix_node = parent_node;
			}
			active.last_char_index++;
			active.canonize();
		}
	}
	
	@Override
	public String toString() {
		return "SuffixTree [treeHeight=" + treeHeight +  " nodes=" + Count + "]";
	}

	/**
	 * String must not contain termination character
	 * @param str
	 * @return
	 */
	public int findStringNode(String str){
		byte[] bytes = null;
		if(str.indexOf(endChar)>-1) {
			System.out.println("String " + str + " contains the end marker");
			return -1;
		}
		bytes = convertStringToByteIndex(str);
		return findNodeForPath(bytes);
	}
	
	int findCountForPath(byte[] bytes) {
		int node = findNodeForPath(bytes) ;
		if(node == -1){return 0;}
		return leafCount[node];
	}
	
	int findNodeForPath(byte[] bytes) {
		int nodeIndex = 0;
		Node node = nodes[nodeIndex];
		Edge next = null;
		for (int strIndex = 0; strIndex < T.length; ) {
			byte idx = bytes[strIndex];
			if(idx==alphabet.length){
				return -1;
			}
			next = node.find(idx);
			if(next==null){
				return -1;
			}
			for (int j = next.first_char_index; j <= next.last_char_index; j++) {
				idx = bytes[strIndex];
				if(idx != T[j]){
					return -1;
				}
				//matched. Go to next char (incr strindx)
				strIndex++;
				if(strIndex == bytes.length){
					//full match
					return next.end_node;
				}
			}
			//done with next
			nodeIndex = next.end_node;
			node = nodes[nodeIndex];
			if(node==null){
				return -1;
			}
		}
		return -1;
	}

	public int findCount(String str){
		int node = findStringNode( str);
		if(node == -1){return 0;}
		return leafCount[node];
	}
	/**
	 * Function to print all contents and details of suffix tree
	 **/
	public void dump_edges(int current_n, int start_node, int level) {
		if(level == 0){
			System.out.println(" Start  End  Suf  First Last  String");
			System.out.println(" Node  Node       char  char        ");
		}
		Node current = nodes[start_node];
		for (Edge edge : current) {
				printEdge(current_n, edge,"");
				if (edge.last_char_index>=N){
					// this edge points to a leaf. No need to make the recursive call,
					// just handle the logic from within this call
					continue;
				} 
				dump_edges(current_n, edge.end_node, level+1);
		}
	}
	
	private char getChar(int l){
		if(T[l]==alphabet.length){return endChar;}
		return alphabet[T[l]];
	}
	private void printEdge(int current_n, Edge edge, String message) {
		if(message.length()>0){System.out.println(message);}
		int suffix_node = edge.last_char_index>=N?-5:nodes[edge.end_node].suffix_node;
		System.out.printf("%5d %5d %3d %5d %6d   ", edge.start_node, edge.end_node, suffix_node,
				edge.first_char_index, edge.last_char_index);
		int top = (current_n > edge.last_char_index)?edge.last_char_index:current_n;
		for (int l = edge.first_char_index; l <= top; l++)
			System.out.print(getChar(l));
		System.out.println();
	}

	/**
	 *  Normally, suffix trees require that the last character in the input string be unique.  
	 *  If you don't do this, your tree will contain suffixes that don't end in leaf nodes.  
	 *  This is often a useful requirement. You can build a tree in this program without 
	 *  meeting this requirement, the code will end the end character for you;
	 **/
	public static void main(String[] args) throws IOException {
		/**
		 * BufferedReader br = new BufferedReader(new
		 * InputStreamReader(System.in)); System.out.println("Enter string\n");
		 * String str = br.readLine();
		 **/
		String str = "ABDXABDBDABXAB";
		SuffixTree.endChar = 'C';
		SuffixTree.alphabet = new char[]{'A','X', 'B', 'D'};
		SuffixTree st = new SuffixTree(convertStringToByteIndex(str));
		st.buildTree();
		
		st.dump_edges(st.N, 0, 0);
		//st.validate();
		st.walk_tree_short(0, new StringBuilder(), 0, -1);
		String[] test = new String[] {"A","AB", "BXA", "B","ABABXABCA","ABB","ABDBD", "XADC", str};
		for (String s : test) {
			System.out.println(s + " " + st.findStringNode(s ) + ": count " + st.findCount(s));
		}
	}

	 void buildTree() {
		for (int i = 0; i <= N; i++) {
			addPrefix(active, i);
		}
		if(T[N]!=alphabet.length){
			addPrefix(active, N+1);
		}
		leafCount = new int[Count];
		makeTreeExplicit(0,0);
	}

	static byte[] convertStringToByteIndex(String str) {
		byte[] T =  new byte[str.length()];
		for (int i = 0; i < str.length(); i++) {
			byte idx = 0;
			char c = str.charAt(i);
			if(c == endChar){
				if(i!=str.length()-1){
					throw new IllegalStateException("end char " + c  + " is not at end, " + str);
				}
				idx = (byte) alphabet.length;
			} else {
				for(;idx <alphabet.length;idx++){
					if(alphabet[idx] == c){break;}
				}
				if(idx==alphabet.length){
					throw new IllegalStateException("alphabet  doesn't contain char " + c  + " " + str);
				}
			}
			T[i] = idx;
		}
		return T;
	}
	
	/**
	 * The validation code consists of two routines. All it does is traverse the
	 * entire tree. walk_tree() calls itself recursively, building suffix
	 * strings up as it goes. When walk_tree() reaches a leaf node, it checks to
	 * see if the suffix derived from the tree matches the suffix starting at
	 * the same point in the input text. If so, it tags that suffix as correct
	 * in the GoodSuffixes[] array. When the tree has been traversed, every
	 * entry in the GoodSuffixes array should have a value of 1.
	 * 
	 * In addition, the BranchCount[] array is updated while the tree is walked
	 * as well. Every count in the array has the number of child edges emanating
	 * from that node. If the node is a leaf node, the value is set to -1. When
	 * the routine finishes, every node should be a branch or a leaf. The number
	 * of leaf nodes should match the number of suffixes (the length) of the
	 * input string. The total number of branches from all nodes should match
	 * the node count.
	 **/
	int[] GoodSuffixes;
	int[] BranchCount;

	void validate() {
		GoodSuffixes = new int[N+1];
		BranchCount = new int[Count];
		walk_tree(0, 0, new StringBuilder(), 0, -1);
		int error = 0;
		for (int i = 0; i < N; i++){
			if (GoodSuffixes[i] != 1) {
				System.out.println("Suffix " + i + " count wrong! " + GoodSuffixes[i]);
				error++;
			}
		}
		if (error == 0)
			System.out.println("All Suffixes present!");
		int leaf_count = 0;
		int branch_count = 0;
		for (int i = 0; i < Count; i++) {
			if (BranchCount[i] == 0)
				System.out.println("Logic error on node " + i + ", not a leaf or internal node!\n");
			else if (BranchCount[i] == -1)
				leaf_count++;
			else
				branch_count += BranchCount[i];
			if(nodes[i]!=null && nodes[i].implicitLeafCount>0){
				leaf_count++;
			}
		}
		int leaves = T[N]!=alphabet.length?N+2:N+1;
		System.out.println("Leaf count : " + leaf_count + (leaf_count == (leaves) ? " OK" : " Error!"));
		System.out.println("Node 0 Leaf count : " + leafCount[0] + (leafCount[0] == (leaves) ? " OK" : " Error!"));
		System.out.println("Branch count : " + branch_count + (branch_count == (Count - 1) ? " OK" : " Error!"));
	}

	void walk_tree(int start_node, int last_char_so_far, StringBuilder CurrentString, int level, int maxLevel) {
		printTreeNode(start_node,  CurrentString, level, maxLevel, false);
		Node current = nodes[start_node];
		for (Edge edge : current) {
			if(GoodSuffixes!=null){
				if(nodes[start_node].implicitLeafCount>0 && GoodSuffixes[last_char_so_far ]==0){
					GoodSuffixes[last_char_so_far ]++;
				}
				if (BranchCount[edge.start_node] < 0) {
					System.out.println("in non-leaf section, Logic error on node " + edge.start_node);
				}
				BranchCount[edge.start_node]++;
			}
			int l = last_char_so_far;
			
			for (int j = edge.first_char_index; j <= edge.last_char_index; j++) {
				CurrentString.insert(l, getChar(j));
				l++;
			}

			CurrentString = new StringBuilder(CurrentString.substring(0, l));

			if (edge.last_char_index>=N && ( (maxLevel==-1) || level<maxLevel)){
				// this edge points to a leaf. No need to make the recursive call,
				// just handle the logic from within this call
				printTreeNode(edge.end_node,  CurrentString, level+1, maxLevel, false);
				if(GoodSuffixes!=null){
					if (BranchCount[edge.end_node] > 0) {
						System.out.println("Logic error on (leaf) node " + edge.end_node);
					}
					BranchCount[edge.end_node]--;
					int idx = T[N]!=alphabet.length?CurrentString.length():CurrentString.length() - 1;
					if(idx<=N){GoodSuffixes[idx]++;}
				}
				continue;
			} 
			if( (maxLevel==-1) || level<maxLevel){
				walk_tree(edge.end_node, l, CurrentString, level+1, maxLevel);
			}
		}
		if(level==0){
			System.out.println();
		}
	}

	void walk_tree_short(int start_node,  StringBuilder CurrentString, 
			int level, int maxLevel ) {
		walk_tree_short( start_node,   CurrentString, 
				 level,  maxLevel, false );
	}
	
	void walk_tree_short(int start_node,  StringBuilder CurrentString, 
			int level, int maxLevel, boolean lastChild ) {
		printTreeNode(start_node,  CurrentString, level, maxLevel, lastChild);
		Node current = nodes[start_node];
		Iterator<Edge> iterator = current.iterator();
		while (iterator.hasNext()) {
			Edge edge = iterator.next();
			CurrentString = new StringBuilder();
			int count = 0;
			for (int j = edge.first_char_index; j <= edge.last_char_index; j++) {
				CurrentString.append( getChar(j));
				++count;
				if((maxLevel>-1) && count>maxLevel && j < edge.last_char_index){
					CurrentString.append("..");
					break;
				}
			}

			if (edge.last_char_index>=N && ( (maxLevel==-1) || level<maxLevel)){
				// this edge points to a leaf. No need to make the recursive call,
				// just handle the logic from within this call
				printTreeNode(edge.end_node, CurrentString, level+1, maxLevel,!iterator.hasNext());
				continue;
			} 
			if( (maxLevel==-1) || level<maxLevel){
				walk_tree_short(edge.end_node,  CurrentString, level+1, maxLevel,!iterator.hasNext());
			}
		}
		if(level==0){
			System.out.println();
		}
	}

	private void printTreeNode(int start_node,  StringBuilder CurrentString, 
			int level, int maxLevel, boolean lastChild) {
		System.out.println();
		for (int i = 0; i < level; i++) {
			if(i==level-1){
				   System.out.print("***" );
			} else if(lastChild && i==level-2){
				   System.out.print(">>>" );
			} else {
				   System.out.print("  |" );
			}
		}
		String out = "";
		if(maxLevel>-1){		
			int max = ((CurrentString.length())<=(maxLevel+1))? CurrentString.length():(maxLevel+1);
			out = CurrentString.substring(0, max);
			if(CurrentString.length()>max){out += "..";}
		}else {
			out = CurrentString.toString();
		}
		System.out.print(" <"+ intFormat.format(start_node) +  "> (" + leafCount[start_node] +  ") [" + out + "]" );
	}
}