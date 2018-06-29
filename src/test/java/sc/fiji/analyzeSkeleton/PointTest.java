
package sc.fiji.analyzeSkeleton;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import org.junit.Test;

/**
 * Tests for {@link Point}
 *
 * @author Richard Domander (Royal Veterinary College, London)
 */
public class PointTest {

	@Test
	public void testClone() throws Exception {
		final Point point = new Point(1, 2, 3);

		final Point clone = point.clone();

		assertEquals(point, clone);
		assertFalse(point == clone);
	}

}
