package comp559.lcp;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;

import java.util.Objects;

import javax.vecmath.GMatrix;
import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;

/**
 * Implementation of a contact constraint.
 * @author kry
 */
public class Contact {

    /** Next available contact index, used for determining which rows of the jacobian a contact uses */
    static public int nextContactIndex = 0;
    
    /** Index of this contact, determines its rows in the jacobian */
    int index;
    
    /** First RigidBody in contact */
    RigidBody body1;
    
    /** Second RigidBody in contact */
    RigidBody body2;
    
    /** Contact normal in world coordinates */
    Vector2d normal = new Vector2d();
    
    /** Position of contact point in world coordinates */
    Point2d contactW = new Point2d();


    double[] jacobianRowOne;
    double[] jacobianRowTwo;
    double[] massMat;
    
    /**
     * Creates a new contact, and assigns it an index
     * @param body1
     * @param body2
     * @param contactW
     * @param normal
     */
    public Contact( RigidBody body1, RigidBody body2, Point2d contactW, Vector2d normal ) {
        this.body1 = body1;
        this.body2 = body2;
        this.contactW.set( contactW );
        this.normal.set( normal );        
        index = nextContactIndex++;        
        // TODO: Objective 3 - Compute and store the contact Jacobian
        
        Vector2d r1 = new Vector2d();
        Vector2d r2 = new Vector2d();
        r1.sub(contactW, body1.x);
        r2.sub(contactW, body2.x);
        jacobianRowOne = new double[6];
        jacobianRowOne[0] = -normal.x;
        jacobianRowOne[1] = -normal.y;
        jacobianRowOne[2] = -r1.x*normal.y + r1.y*normal.x;
        jacobianRowOne[3] = normal.x;
        jacobianRowOne[4] = normal.y;
        jacobianRowOne[5] = r2.x*normal.y - r2.y*normal.x;
        
        Vector2d tan = new Vector2d(-normal.y, normal.x);
        jacobianRowTwo = new double[6];
        jacobianRowTwo[0] = -tan.x;
        jacobianRowTwo[1] = -tan.y;
        jacobianRowTwo[2] = -r1.x*tan.y + r1.y*tan.x;
        jacobianRowTwo[3] = tan.x;
        jacobianRowTwo[4] = tan.y;
        jacobianRowTwo[5] = r2.x*tan.y - r2.y*tan.x;
        
        massMat = new double[] {body1.minv, body1.minv, body1.jinv, 
        		body2.minv, body2.minv, body2.jinv};
    }
    
    /**
     * Draws the contact points
     * @param drawable
     */
    public void display( GLAutoDrawable drawable ) {
        GL2 gl = drawable.getGL().getGL2();
        gl.glPointSize(3);
        gl.glColor3f(.7f,0,0);
        gl.glBegin( GL.GL_POINTS );
        gl.glVertex2d(contactW.x, contactW.y);
        gl.glEnd();
    }
    
    /**
     * Draws the connections between bodies to visualize the 
     * the adjacency structure of the matrix as a graph.
     * @param drawable
     */
    public void displayConnection( GLAutoDrawable drawable ) {
        GL2 gl = drawable.getGL().getGL2();
        // draw a line between the two bodies but only if they're both not pinned
        if ( !body1.pinned && ! body2.pinned ) {
            gl.glLineWidth(2);
            gl.glColor4f(0,.3f,0, 0.5f);
            gl.glBegin( GL.GL_LINES );
            gl.glVertex2d(body1.x.x, body1.x.y);
            gl.glVertex2d(body2.x.x, body2.x.y);
            gl.glEnd();
        }
    }
    
    
    @Override
    public boolean equals(Object o) {
	  if (this == o) return true;
	  if (o == null || getClass() != o.getClass()) return false;
	  Contact contact = (Contact) o;
	  return (body1.equals(contact.body1) && body2.equals(contact.body2)) || (body2.equals(contact.body1) && body1.equals(contact.body1));
    }

    @Override
    public int hashCode() {
    	return body1.hashCode() + body2.hashCode();
//      return Objects.hash(body1, body2);
    }
    
}
