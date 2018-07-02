package sc.fiji.analyzeSkeleton;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;

import org.junit.Test;

/**
 * Tests for {@link Edge}
 *
 * @author Richard Domander (Royal Veterinary College, London)
 */
public class EdgeTest {

    @Test
    public void testCloneUnconnected() throws Exception {
        final ArrayList<Point> slabs = new ArrayList<>();
        slabs.add(new Point(1, 2, 3));
        slabs.add(new Point(4, 5, 6));
        final Edge edge = new Edge(null, null, slabs, 7.0, 8.0, 9.0, 10.0);
        edge.setType(11);

        final Edge clone = edge.clone(null, null);

        assertTrue(edge != clone);
        assertEquals(edge.getColor(), clone.getColor(), 1e-12);
        assertEquals(edge.getColor3rd(), clone.getColor3rd(), 1e-12);
        assertEquals(edge.getLength(), clone.getLength(), 1e-12);
        assertEquals(edge.getLength_ra(), clone.getLength_ra(), 1e-12);
        assertEquals(edge.getType(), clone.getType());
        final ArrayList<Point> cloneSlabs = clone.getSlabs();
        assertEquals(slabs.size(), cloneSlabs.size());
        for (int i = 0; i < slabs.size(); i++) {
            assertEquals(slabs.get(i), cloneSlabs.get(i));
        }
    }

}