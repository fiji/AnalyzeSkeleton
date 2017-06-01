
package sc.fiji.analyzeSkeleton;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;

import org.junit.Test;

/**
 * Tests for {@link Vertex}
 *
 * @author Richard Domander
 */
public class VertexTest {

	@Test
	public void testCloneUnconnected() throws Exception {
		final Vertex vertex = new Vertex();
		vertex.setVisited(true, 4);
		vertex.addPoint(new Point(1, 2, 3));
		vertex.addPoint(new Point(4, 5, 6));

		final Vertex clone = vertex.cloneUnconnected();

		assertTrue(clone != vertex);
		assertEquals(vertex.isVisited(), clone.isVisited());
		assertEquals(vertex.getVisitOrder(), clone.getVisitOrder());
		final ArrayList<Point> points = vertex.getPoints();
		final ArrayList<Point> clonePoints = clone.getPoints();
		assertEquals(points.size(), clonePoints.size());
		for (int i = 0; i < points.size(); i++) {
			assertEquals(points.get(i), clonePoints.get(i));
		}
	}

}
