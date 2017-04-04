/**
 * 
 */
package suffix;

import static org.junit.Assert.*;

import org.junit.Test;

/**
 * @author oshanker
 *
 */
public class SuffixTreeTest {

	/**
	 * Test method for {@link suffix.SuffixTree#buildTree()}.
	 */
	@Test
	public void testBuildTree() {
		byte[] T = new byte[1000];
		SuffixTree.endChar = 'Z';
		SuffixTree.alphabet = new char[]{'A', 'B', 'C','D'};
		for (int i = 0; i < T.length; i++) {
			T[i] = (byte) (Math.random()*SuffixTree.alphabet.length);
		}
		SuffixTree st = new SuffixTree(T);
		st.buildTree();
		
		st.walk_tree_short(0, new StringBuilder(), 0, 3);
	}

}
