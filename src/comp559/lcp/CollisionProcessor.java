package comp559.lcp;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;

import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.CollapsiblePanel;
import mintools.swing.VerticalFlowPanel;

/**
 * Class for detecting and resolving collisions.  Currently this class uses penalty forces between rigid bodies.
 * @author kry
 */
public class CollisionProcessor {

    private List<RigidBody> bodies;
    
    /**
     * The current contacts that resulted in the last call to process collisions
     */
    public ArrayList<Contact> contacts = new ArrayList<Contact>();
    
    /**
     * Creates this collision processor with the provided set of bodies
     * @param bodies
     */
    public CollisionProcessor( List<RigidBody> bodies ) {
        this.bodies = bodies;
    }
    
    /** keeps track of the time used for collision detection on the last call */
    double collisionDetectTime = 0;
    
    /** keeps track of the time used to solve the LCP based velocity update on the last call */
    double collisionSolveTime = 0;
    
    public Map<Contact, double[]> pair = new HashMap<>();
    
    public void createU(double[] u, Contact ct) {
    	u[0] = ct.body1.v.x;
        u[1] = ct.body1.v.y;
        u[2] = ct.body1.omega;
        u[3] = ct.body2.v.x;
        u[4] = ct.body2.v.y;
        u[5] = ct.body2.omega;
    }
    
    public void createF(double[] f, Contact ct) {
    	f[0] = ct.body1.force.x;
        f[1] = ct.body1.force.y;
        f[2] = ct.body1.torque;
        f[3] = ct.body2.force.x;
        f[4] = ct.body2.force.y;
        f[5] = ct.body2.torque;
    }
    
    public double[] getconstraints(double dt, Contact ct, double[] u, double[] f) {
    	double Nconstraint = 0;
        double Tconstraint = 0;
        double Rconstraint = 0;
        
    	for(int i = 0; i < 6; i++) {
    		Nconstraint += ct.jacobianRowOne[i] * (u[i] + dt * ct.massMat[i] * f[i]);
    		Tconstraint += ct.jacobianRowTwo[i] * (u[i] + dt * ct.massMat[i] * f[i]);
    		Rconstraint += ct.jacobianRowOne[i] * restitution.getValue() * u[i];
    	}
    	return new double[] {Nconstraint, Tconstraint, Rconstraint};
    }
    
    public void assembleB(double dt, double[] b) {
    	for(Contact ct: contacts) {
            double[] u = new double[6];
            createU(u, ct);
            double[] f = new double[6];
            createF(f, ct);
            double[] constraints = getconstraints(dt, ct, u, f);

        	b[2 * ct.index] = constraints[0] + constraints[2];
        	b[2 * ct.index + 1] = constraints[1];
        }
    }
    
    public void assembleD(double dt, double[] D) {
    	for(Contact ct: contacts) {
        	double d1 = 0;
        	double d2 = 0;
        	for(int i = 0; i < 6; i++) {
        		d1 += ct.jacobianRowOne[i] * ct.jacobianRowOne[i] * ct.massMat[i];
        		d2 += ct.jacobianRowTwo[i] * ct.jacobianRowTwo[i] * ct.massMat[i];
        	}
        	D[2 * ct.index] = d1;
        	D[2 * ct.index + 1] = d2;
        }
    }
    
    public void normalizeB(double[] b, double[] D, double[] bprime) {
    	for(int i = 0; i < 2*contacts.size(); i++) {
        	bprime[i] = b[i]/D[i];
        }
    }
    
    public void setLambdaI(double[] lambdaI, Contact ctct, double[] tmp) {
    	lambdaI[2 * ctct.index] = tmp[6];
        lambdaI[2 * ctct.index + 1] = tmp[7];
    }
    
    public void settmpdeltav1(double[] T, Contact ctct) {
    	for(int i = 0; i < 6; i++) {
			T[i] = ctct.massMat[i] * ctct.jacobianRowOne[i];
		}
    }
    
    public void settempdeltav2(double[] T, Contact ctct) {
    	for(int i = 0; i < 6; i++) {
			T[i] = ctct.massMat[i] * ctct.jacobianRowTwo[i];
		}
    }
    
    /**
     * Processes all collisions 
     * @param dt time step
     */
    public void processCollisions( double dt ) {
        contacts.clear();
        Contact.nextContactIndex = 0;
        
        long now = System.nanoTime();
        broadPhase();
        collisionDetectTime = ( System.nanoTime() - now ) * 1e-9;
                
        if ( contacts.size() > 0  && doLCP.getValue() ) {
            now = System.nanoTime();
            double mu = friction.getValue();

            // TODO: Objective 3 - Compute velocity update with iterative solve of contact constraint matrix.
            
            int iterN = iterations.getValue();
            int n = contacts.size();
            int N = bodies.size();
            
            //construct D and b:
            double[] D = new double[2*n];
            double[] b = new double[2*n];
            
            assembleB(dt, b);
            assembleD(dt, D);
            
            double[] bprime = new double[2*n];
            normalizeB(b, D, bprime);
            
            //compute lambda recursively (projective Gauss Seidel)
            double[] lambdaI = new double[2*n];
            double[] deltaV = new double[3*N];

            if (warmStart.getValue()) {
	            for(Contact ctct: contacts) {
	            	double[] tmp = pair.remove(ctct);
	            	if(tmp!=null) {
                    	for(int i = 0; i < 3; i++) {
                			deltaV[ctct.body1.index * 3 + i] = tmp[i];
                		}
                    	for(int i = 0; i < 3; i++) {
                    		deltaV[ctct.body2.index * 3 + i] = tmp[i+3];
                    	}
                        setLambdaI(lambdaI, ctct, tmp);
	            	}
	            }
            }
            
            for(int j = 0; j < iterN; j++) {
            	if(useShuffle.getValue()) Collections.shuffle(contacts);
            	for(Contact ctct: contacts) {
            		int idx = 2*ctct.index;
            		int body1_idx = ctct.body1.index * 3;
            		int body2_idx = ctct.body2.index * 3;
            		
            		double lambdaN = lambdaI[idx] - bprime[idx];
            		double lambdaF = lambdaI[idx+1] - bprime[idx+1];
            		for(int i = 0; i < 3; i++) {
            			lambdaN -= ctct.jacobianRowOne[i] / D[idx] * deltaV[body1_idx + i];
            			lambdaN -= ctct.jacobianRowOne[i+3] / D[idx] * deltaV[body2_idx + i];
            			
            			lambdaF -= ctct.jacobianRowTwo[i] / D[idx+1] * deltaV[body1_idx + i];
            			lambdaF -= ctct.jacobianRowTwo[i+3] / D[idx+1] * deltaV[body2_idx + i];
            		}
            		lambdaN = Math.max(lambdaN, 0);
            		double deltalambda1 = lambdaN - lambdaI[idx];
            		lambdaI[idx] = lambdaN; //update
            		
            		double[] T = new double[6];
            		settmpdeltav1(T, ctct);
            		
            		for(int i = 0; i < 3; i++) {
            			deltaV[body1_idx + i] += T[i]*deltalambda1;
            			deltaV[body2_idx + i] += T[i+3]*deltalambda1;
            		}
            		
            		
            		double lambdaT = mu * lambdaN;
            		lambdaF = Math.min(lambdaT, Math.max(-lambdaT, lambdaF));  //projection
            		double deltalambda2 = lambdaF - lambdaI[idx+1];
            		lambdaI[idx+1] = lambdaF;
            		settempdeltav2(T, ctct);
            		for(int i = 0; i < 3; i++) {
            			deltaV[body1_idx + i] += T[i]*deltalambda2;
            			deltaV[body2_idx + i] += T[i+3]*deltalambda2;
            		}          		
            	}
            	
            }
            
            if (warmStart.getValue()){
            	pair.clear();
            	for(Contact ctct: contacts) {
            		int idx = 2*ctct.index;
            		double lambda1 = lambdaI[idx];
                	double lambda2 = lambdaI[idx + 1];
                	
                	double[] tempdeltav = new double[6];
                	for(int i = 0; i < 6; i++) {
                		tempdeltav[i] = ctct.massMat[i] * (ctct.jacobianRowOne[i]*lambda1 + ctct.jacobianRowTwo[i]*lambda2);
                	}
                    pair.put(ctct, new double[] {tempdeltav[0], tempdeltav[1], tempdeltav[2], tempdeltav[3], tempdeltav[4], tempdeltav[5], lambda1, lambda2 });
            	}
    		}
            
            for(RigidBody body: bodies) {
            	body.v.x += deltaV[body.index*3];
            	body.v.y += deltaV[body.index*3 + 1];
            	body.omega += deltaV[body.index*3 + 2];
            }
                        
            collisionSolveTime = (System.nanoTime() - now) * 1e-9;
        }
    }
    
    /**
     * Checks for collisions between bodies.  Note that you can optionaly implement some broad
     * phase test such as spatial hashing to reduce the n squared body-body tests.
     * Currently this does the naive n squared collision check.
     */
    private void broadPhase() {
        // Naive n squared body test.. might not be that bad for small number of bodies 
        visitID++;
        for ( RigidBody b1 : bodies ) {
            for ( RigidBody b2 : bodies ) { // not so inefficient given the continue on the next line
                if ( b1.index >= b2.index ) continue;
                if ( b1.pinned && b2.pinned ) continue;                
                narrowPhase( b1, b2 );                
            }
        }        
    }
    
    /**
     * Checks for collision between boundary blocks on two rigid bodies.
     * TODO: Objective 2 - This needs to be improved as the n-squared block test is too slow!
     * @param body1
     * @param body2
     */
    private void narrowPhase( RigidBody body1, RigidBody body2 ) {
        if ( ! useBVTree.getValue() ) {
            for ( Block b1 : body1.blocks ) {
                for ( Block b2 : body2.blocks ) {
                    processCollision( body1, b1, body2, b2 );
                }
            }
        } else {
            // TODO: Objective 2 - Implement code to use hierarchical collision detection on body pairs
            BVNode n1 = body1.root;
            BVNode n2 = body2.root;
            narrow_hierarchy(n1, n2, body1, body2);

        }
    }
    
    public void narrow_hierarchy(BVNode n1, BVNode n2, RigidBody body1, RigidBody body2) {
    	n1.boundingDisc.updatecW();
        n2.boundingDisc.updatecW();
        if(n1.boundingDisc.intersects(n2.boundingDisc)) {
        	traverse(body1, n1, body2, n2);
        }
    }

    
    
    private void trytraverse1(RigidBody rb1, BVNode n1, RigidBody rb2, BVNode n2) {
    	if (n1.boundingDisc.intersects(n2.boundingDisc)) {
            traverse(rb1, n1, rb2, n2);
        }
    }
    
    private void trytraverse2(RigidBody rb1, BVNode n1, RigidBody rb2, BVNode n2) {
    	if (n2.boundingDisc.intersects(n1.boundingDisc)) {
            traverse(rb1, n1, rb2, n2);
        }
    }
    
    private void update_LR_visit_ID(BVNode L, BVNode R, int ID) {
    	if (L.visitID != ID) {
        	L.visitID = ID;
            L.boundingDisc.updatecW();
        }
    	if (R.visitID != ID) {
        	R.visitID = ID;
            R.boundingDisc.updatecW();
        }
    }
    
    private boolean bothLeaf(BVNode Bvn_1, BVNode Bvn_2) {
    	return Bvn_1.isLeaf() && Bvn_2.isLeaf();
    }

    private void traverse(RigidBody rb1, BVNode Bvn_1, RigidBody rb2, BVNode Bvn_2){
        if (bothLeaf(Bvn_1, Bvn_2)) {
            processCollision(rb1, Bvn_1.leafBlock, rb2, Bvn_2.leafBlock);
        } 
        else if (Bvn_1.isLeaf() && !Bvn_2.isLeaf()) {
            update_LR_visit_ID(Bvn_2.child1, Bvn_2.child2, visitID);
            trytraverse1(rb1, Bvn_1, rb2, Bvn_2.child1);
            trytraverse1(rb1, Bvn_1, rb2, Bvn_2.child2);
        } 
        else if (!Bvn_1.isLeaf() && Bvn_2.isLeaf()) {
            update_LR_visit_ID(Bvn_1.child1, Bvn_1.child2, visitID);
            trytraverse2(rb1, Bvn_1.child1, rb2, Bvn_2);
            trytraverse2(rb1, Bvn_1.child2, rb2, Bvn_2);
        } 
        else if (Bvn_1.boundingDisc.r > Bvn_2.boundingDisc.r) {
        	update_LR_visit_ID(Bvn_1.child1, Bvn_1.child2, visitID);
            trytraverse2(rb1, Bvn_1.child1, rb2, Bvn_2);
            trytraverse2(rb1, Bvn_1.child2, rb2, Bvn_2);
        }
        else {
            update_LR_visit_ID(Bvn_2.child1, Bvn_2.child2, visitID);
            trytraverse1(rb1, Bvn_1, rb2, Bvn_2.child1);
            trytraverse1(rb1, Bvn_1, rb2, Bvn_2.child2);
        }

    }
    
    

    /** 
     * The visitID is used to tag boundary volumes that are visited in 
     * a given time step.  Marking boundary volume nodes as visited during
     * a time step allows for a visualization of those used, but it can also
     * be used to more efficiently update the centeres of bounding volumes
     * (i.e., call a BVNode's updatecW method at most once on any given timestep)
     */
    int visitID = 0;
    
    /**
     * Resets the state of the collision processor by clearing all
     * currently identified contacts, and reseting the visitID for
     * tracking the bounding volumes used
     */
    public void reset() {
        contacts.clear();
        Contact.nextContactIndex = 0;
        visitID = 0;            
    }
    
    // some working variables for processing collisions
    private Point2d tmp1 = new Point2d();
    private Point2d tmp2 = new Point2d();
    private Point2d contactW = new Point2d();
    private Vector2d force = new Vector2d();
    private Vector2d contactV1 = new Vector2d();
    private Vector2d contactV2 = new Vector2d();
    private Vector2d relativeVelocity = new Vector2d();
    private Vector2d normal = new Vector2d();
        
    /**
     * Processes a collision between two bodies for two given blocks that are colliding.
     * Currently this implements a penalty force
     * @param body1
     * @param b1
     * @param body2
     * @param b2
     */
    private void processCollision( RigidBody body1, Block b1, RigidBody body2, Block b2 ) {        
        double k = contactSpringStiffness.getValue();
        double c1 = contactSpringDamping.getValue();
        double threshold = separationVelocityThreshold.getValue();
        boolean useSpring = enableContactSpring.getValue();
        boolean useDamping = enableContactDamping.getValue();
        
        body1.transformB2W.transform( b1.pB, tmp1 );
        body2.transformB2W.transform( b2.pB, tmp2 );
        double distance = tmp1.distance(tmp2);
        if ( distance < Block.radius * 2 ) {
            // contact point at halfway between points 
            // NOTE: this assumes that the two blocks have the same radius!
            contactW.interpolate( tmp1, tmp2, .5 );
            // contact normal
            normal.sub( tmp2, tmp1 );
            normal.normalize();
            // create the contact
            Contact contact = new Contact( body1, body2, contactW, normal);
            // simple option... add to contact list...
            
            contacts.add( contact );
            
            
            
            if ( ! doLCP.getValue() ) {
                // compute relative body velocity at contact point
                body1.getSpatialVelocity( contactW, contactV1 );
                body2.getSpatialVelocity( contactW, contactV2 );
                relativeVelocity.sub( contactV1, contactV2 );
                if ( -relativeVelocity.dot( normal ) < threshold ) {
                    if ( useSpring ) {
                        // spring force
                        double interpenetration = distance - Block.radius * 2; // a negative quantity
                        force.scale( -interpenetration * k, normal );
                        body2.applyContactForceW(contactW, force);
                        force.scale(-1);
                        body1.applyContactForceW(contactW, force);
                    }
                    if ( useDamping ) {
                        // spring damping forces!
                        // vertical
                        force.scale( relativeVelocity.dot(normal) * c1, normal );                    
                        body2.applyContactForceW( contactW, force );
                        force.scale(-1);
                        body1.applyContactForceW( contactW, force );
                    }
                }
            }
        }
    }
   
    /** Stiffness of the contact penalty spring */
    private DoubleParameter contactSpringStiffness = new DoubleParameter("penalty contact stiffness", 1e3, 1, 1e5 );
    
    /** Viscous damping coefficient for the contact penalty spring */
    private DoubleParameter contactSpringDamping = new DoubleParameter("penalty contact damping", 10, 1, 1e4 );
    
    /** Threshold for the relative velocity in the normal direction, for determining if spring force will be applied. */
    private DoubleParameter separationVelocityThreshold = new DoubleParameter( "penalty separation velocity threshold (controls bounce)", 1e-9, 1e-9, 1e3 );
    
    /** Enables the contact penalty spring */
    private BooleanParameter enableContactSpring = new BooleanParameter("enable penalty contact spring", true );
    
    /** Enables damping of the contact penalty spring */
    private BooleanParameter enableContactDamping = new BooleanParameter("enable penalty contact damping", true );
    
    /** Restitution parameter for contact constraints */
    public DoubleParameter restitution = new DoubleParameter( "restitution (bounce)", 0, 0, 1 );
    
    /** Coulomb friction coefficient for contact constraint */
    public DoubleParameter friction = new DoubleParameter("Coulomb friction", 0.33, 0, 2 );
    
    /** Number of iterations to use in projected Gauss Seidel solve */
    public IntParameter iterations = new IntParameter("iterations for GS solve", 10, 1, 500);
    
    /** Flag for switching between penalty based contact and contact constraints */
    private BooleanParameter doLCP = new BooleanParameter( "do LCP solve", false );
    
    /** Flag for enabling the use of hierarchical collision detection for body pairs */
    private BooleanParameter useBVTree = new BooleanParameter( "use BVTree", false );

    /** Flag for enabling randomization */
    private BooleanParameter useShuffle = new BooleanParameter("use Shuffle", false);

    private BooleanParameter warmStart = new BooleanParameter("use Warm Start", false);

    /**
     * @return controls for the collision processor
     */
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        vfp.setBorder( new TitledBorder("Collision Processing Controls") );
        vfp.add( useBVTree.getControls() );
        vfp.add( useShuffle.getControls());
        vfp.add( doLCP.getControls() );
        vfp.add( warmStart.getControls());
        vfp.add( iterations.getSliderControls() );
        vfp.add( restitution.getSliderControls(false) );
        vfp.add( friction.getSliderControls(false) );
        
        VerticalFlowPanel vfp2 = new VerticalFlowPanel();
        vfp2.setBorder( new TitledBorder("penalty method controls") );
        vfp2.add( contactSpringStiffness.getSliderControls(true) );
        vfp2.add( contactSpringDamping.getSliderControls(true) );
        vfp2.add( separationVelocityThreshold.getSliderControls( true ) );
        vfp2.add( enableContactDamping.getControls() );
        vfp2.add( enableContactSpring.getControls() );
        
        CollapsiblePanel cp = new CollapsiblePanel(vfp2.getPanel());
        cp.collapse();
        vfp.add( cp );        
        return vfp.getPanel();
    }
    
}
