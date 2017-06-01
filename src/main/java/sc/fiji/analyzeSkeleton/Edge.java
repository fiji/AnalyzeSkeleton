/*
 * #%L
 * AnalyzeSkeleton_ plugin for ImageJ.
 * %%
 * Copyright (C) 2008 - 2017 Ignacio Arganda-Carreras.
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package sc.fiji.analyzeSkeleton;

import java.util.ArrayList;
import java.util.stream.Collectors;

/**
 * This class represents the edge between two vertices in an undirected graph.
 */
public class Edge
{
	/** "tree" edge classification constant for Depth-first search (DFS) */
	public final static int TREE = 0;
	/** "back" edge classification constant for Depth-first search (DFS) */
	public final static int BACK = 1;
	/** not yet defined edge classification constant for Depth-first search (DFS) */
	public final static int UNDEFINED = -1;

	/** DFS classification */
	private int type = Edge.UNDEFINED;

	/** vertex at one extreme of the edge */
	private Vertex v1 = null;
	/** vertex at the other extreme of the edge */
	private Vertex v2 = null;
	/** list of slab voxels belonging to this edge */
	private ArrayList <Point> slabs = null;
	/** length of the edge */
	private double length = 0;
	/** average color of edge */
	private double color = 0;
	/** average color of inner third of  edge */
	private double color3rd = 0;
	/** length calculated by running average ovr 5 Pixel */
	private double length_ra = 0;


	/**
	 * Create an edge of specific vertices and list of slab voxels.
	 * @param v1 first vertex
	 * @param v2 second vertex
	 * @param slabs list of slab voxels
	 * @param length calibrated edge length
	 */
	public Edge(
			Vertex v1, 
			Vertex v2, 
			ArrayList<Point> slabs,
			double length)
	{
		this.v1 = v1;
		this.v2 = v2;
		this.slabs = slabs;
		this.length = length;
	}

	/**
	 * Create an edge of specific vertices and list of slab voxels.
	 * @param v1 first vertex
	 * @param v2 second vertex
	 * @param slabs list of slab voxels
	 * @param length calibrated edge length
	 * @param color3rd average color value of the inner third
	 * @param color average color value
	 * @param length_ra calibrated edge length calculated with running average ofer 5 Pixel
	 */
	public Edge(
			Vertex v1, 
			Vertex v2, 
			ArrayList<Point> slabs,
			double length,
			double color3rd,
			double color,
			double length_ra)
	{
		this.v1 = v1;
		this.v2 = v2;
		this.slabs = slabs;
		this.length = length;
		this.color = color;
		this.color3rd = color3rd;
		this.length_ra = length_ra;

	}

	/**
	 * Get first vertex. 
	 * @return first vertex of the edge
	 */
	public Vertex getV1()
	{
		return this.v1;
	}
	/**
	 * Get second vertex.
	 * @return second vertex of the edge
	 */
	public Vertex getV2()
	{
		return this.v2;
	}
	/**
	 * Get list of slab voxels belonging to the edge.
	 * @return list of slab voxels
	 */
	public ArrayList<Point> getSlabs()
	{
		return this.slabs;
	}
	/**
	 * Set DFS type (BACK or TREE)
	 * @param type DFS classification (BACK or TREE)
	 */
	public void setType(int type)
	{
		this.type = type;
	}
	/**
	 * Get DFS edge type
	 * @return DFS classification type
	 */
	public int getType()
	{
		return this.type;
	}
	/**
	 * Get opposite vertex from a given one.
	 * @param v input vertex
	 * @return opposite vertex in the edge
	 */
	public Vertex getOppositeVertex(Vertex v)
	{
		if(this.v1.equals(v))
			return this.v2;
		else if (this.v2.equals(v))
			return this.v1;
		else 
			return null;
	}

	/**
	 * Set edge length
	 * @param length calibrated edge length
	 */
	public void setLength(double length)
	{
		this.length = length;
	}

	/**
	 * Get edge length
	 * @return calibrated edge length
	 */
	public double getLength()
	{
		return this.length;
	}

	/**
	 * Get edge length_ra (running average)
	 * @return calibrated edge length (running average)
	 */
	public double getLength_ra()
	{
		return this.length_ra;
	}


	public void setColor(double color)
	{
		this.color = color;
	}


	public double getColor()
	{
		return this.color;
	}

	public void setColor3rd(double color)
	{
		this.color3rd = color;
	}


	public double getColor3rd()
	{
		return this.color3rd;
	}

	/**
	 * Clones the Edge with all its properties
	 * <p>
	 * NB Does not clone the vertices!
	 * </p>
	 * 
	 * @param v1 One endpoint of the edge, can be null
	 * @param v2 The other endpoint of the edge, can be null
	 */
	public Edge clone(final Vertex v1, final Vertex v2) {
		final ArrayList<Point> clonedSlabs;
		if (slabs != null) {
			clonedSlabs = slabs.stream().map(Point::clone).collect(Collectors
				.toCollection(ArrayList::new));
		}
		else {
			clonedSlabs = null;
		}
		final Edge clonedEdge = new Edge(v1, v2, clonedSlabs, length, color3rd,
			color, length_ra);
		clonedEdge.setType(type);
		return clonedEdge;
	}

}// end class Edge
