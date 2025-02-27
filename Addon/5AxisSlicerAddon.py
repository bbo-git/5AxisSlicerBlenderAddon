

import bpy
import math
from math import atan2, pi
#workign theta calc
import bmesh
import mathutils 
import os
import subprocess
import json
from bpy.props import StringProperty, FloatProperty, CollectionProperty, PointerProperty


# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTIBILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

bl_info = {
    "name": "4_5th_axis",
    "author": "Jairo Bambang Oetomo",
    "description": "",
    "blender": (2, 80, 0),
    "version": (0, 0, 1),
    "location": "",
    "warning": "",
    "category": "Generic",
}

# -------------------------------------------------------------------
# (2) Slicing Cubes Data (Collection Property)
# -------------------------------------------------------------------
class SlicingCubeItem(bpy.types.PropertyGroup):
    """Holds a reference to one slicing cube object by name"""
    name: bpy.props.StringProperty(default="SlicingCube")
    isProcessed: bpy.props.BoolProperty(default=False, description="Has this slicing cube been processed?")
    
# -------------------------------------------------------------------
# (3) UIList for Slicing Cubes
# -------------------------------------------------------------------
class SLICINGCUBE_UL_items(bpy.types.UIList):
    """List UI for slicing cubes in the add-on panel"""
    def draw_item(self, context, layout, data, item, icon, active_data, active_propname, index=0, flt_flag=0):
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            layout.label(text=item.name, icon='CUBE')
        elif self.layout_type == 'GRID':
            layout.alignment = 'CENTER'
            layout.label(text="", icon='CUBE')

# -------------------------------------------------------------------
# (4) Add/Remove Slicing Cubes
# -------------------------------------------------------------------
class SLICINGCUBE_OT_add(bpy.types.Operator):
    """Add a 50x50x50 slicing cube"""
    bl_idname = "slicingcube.add"
    bl_label = "Add Slicing Cube"

    def execute(self, context):
        scene = context.scene

        # Create the slicing cube
        bpy.ops.mesh.primitive_cube_add(size=50.0, location=(0, 0, 25))
        new_cube = context.active_object
        new_cube.name = self.unique_name("SlicingCube")
        
        scene.cursor.location = (0,0,0)
        bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
        scene.cursor.location = (0,0,scene.build_plate_distance)
        
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='DESELECT')
        
        obj = new_cube
        mesh = obj.data

        # Create a BMesh object from the mesh
        bm = bmesh.from_edit_mesh(mesh)

        # Deselect all faces
        for face in bm.faces:
            face.select = False

        # Find the face(s) with the highest Z-coordinate
        max_z = -float('inf')
        top_faces = []

        for face in bm.faces:
            # Calculate the average Z-coordinate of the face
            face_z = sum(v.co.z for v in face.verts) / len(face.verts)
            
            if face_z > max_z:
                max_z = face_z
                top_faces = [face]
            elif face_z == max_z:
                top_faces.append(face)

        # Select the top face(s)
        for face in top_faces:
            face.select = True

        # Store the index of the top face in a custom property
        if top_faces:
            new_cube["top_face_index"] = top_faces[0].index

        # Update the mesh to reflect the selection
        bmesh.update_edit_mesh(mesh)

        # Extrude the selected face(s)
        bpy.ops.mesh.extrude_region_move(
            TRANSFORM_OT_translate={
                "value": (0, 0, 15),  # Extrude along the Z-axis by 15 units
                "orient_type": 'GLOBAL',  # Use global orientation
            }
        )
        bpy.ops.mesh.merge(type='CENTER')
        bpy.ops.object.mode_set(mode='OBJECT')
        
        #######

        # Track the slicing cube in slicing_cubes_collection
        item = scene.slicing_cubes_collection.add()
        item.name = new_cube.name
        scene.slicing_cubes_index = len(scene.slicing_cubes_collection) - 1

        self.report({'INFO'}, f"Added slicing cube: {new_cube.name}")
        return {'FINISHED'}

    def unique_name(self, base):
        existing = {o.name for o in bpy.data.objects}
        name = base
        i = 1
        while name in existing:
            name = f"{base}.{i:03d}"
            i += 1
        return name

class SLICINGCUBE_OT_remove(bpy.types.Operator):
    """Remove the selected slicing cube from the scene/list"""
    bl_idname = "slicingcube.remove"
    bl_label = "Remove Slicing Cube"

    @classmethod
    def poll(cls, context):
        scene = context.scene
        return (
            len(scene.slicing_cubes_collection) > 0
            and 0 <= scene.slicing_cubes_index < len(scene.slicing_cubes_collection)
        )

    def execute(self, context):
        scene = context.scene
        idx = scene.slicing_cubes_index
        item = scene.slicing_cubes_collection[idx]

        # Delete object if it still exists
        cube_obj = bpy.data.objects.get(item.name)
        if cube_obj:
            bpy.data.objects.remove(cube_obj, do_unlink=True)

        # Remove from collection
        scene.slicing_cubes_collection.remove(idx)
        scene.slicing_cubes_index = min(idx, len(scene.slicing_cubes_collection)-1)

        self.report({'INFO'}, f"Removed slicing cube: {item.name}")
        return {'FINISHED'}

# -------------------------------------------------------------------
# (5) Calculate Theta (Bottom Face => -Z)
# -------------------------------------------------------------------
import bpy
from math import atan2
from mathutils import Vector, Matrix


class SLICINGCUBE_OT_slice(bpy.types.Operator):
    bl_idname = "slicingcube.calculate_theta"
    bl_label = "Slice Cone"
    
    def set_origin_to_a_axis(self, obj):
        # Ensure the object is selected and active
        bpy.context.view_layer.objects.active = obj
        obj.select_set(True)
        
        bpy.ops.object.transform_apply(location=True, rotation=True, scale=True)
        
        # Update cursor location directly
        bpy.context.scene.cursor.location = (0, 0, bpy.context.scene.build_plate_distance)
        
        # Set the object’s origin to the 3D cursor
        bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
        
        # Update and report success
        bpy.context.view_layer.update()
        
        self.report({'INFO'}, f"Origin of {obj.name} set to A xis")
        
        return {'FINISHED'}
           
    def apply_rotations(self, obj, theta, eta):
        """
        Apply the rotations around the Z and X axes to the object.
        """
        # Create rotation matrices
        axis = (0, 0, 1)  # Z-axis
        rotation = self.rotation_matrix(theta, axis)
        obj.matrix_world = rotation @ obj.matrix_world
        
        axis = (1, 0, 0)  # Z-axis
        rotation = self.rotation_matrix(eta, axis)
        obj.matrix_world = rotation @ obj.matrix_world
        
        self.report({'INFO'}, f"{obj.name}rotations applied")
        
    def invert_rotations(self, obj, theta, eta):
        """
        Invert the rotations applied to the object.
        """
        axis = (1, 0, 0)  # Z-axis
        rotation = self.rotation_matrix(-eta, axis)
        obj.matrix_world = rotation @ obj.matrix_world# Create rotation matrices
        
        axis = (0, 0, 1)  # Z-axis
        rotation = self.rotation_matrix(-theta, axis)
        obj.matrix_world = rotation @ obj.matrix_world
        
        self.report({'INFO'}, f"{obj.name}rotations inverted")
    
    def normalize_angle(self, angle):
        return (angle + math.pi) % (2 * math.pi) - math.pi
       
    def calculate_rotations(self, slicing_cube):
        
        print(f"negative: {bpy.context.scene.favor_negative_a}")
        """
        Calculate theta (rotation around Z-axis) and eta (rotation around X-axis)
        based on the bottom face normal of the slicing cube **without modifying it**.
        """
        # Ensure the bottom face exists
        if len(slicing_cube.data.polygons) < 5:
            self.report({'ERROR'}, "Slicing cube does not have enough polygons!")
            return None

        bottom_face = slicing_cube.data.polygons[4]  # Assuming the bottom face is the 4th one
        local_normal = bottom_face.normal  # Normal in local space

        # Transform normal to world space (copy, don't modify the object)
        normal_in_world_space = slicing_cube.matrix_world.to_3x3() @ local_normal
        normal_in_world_space = normal_in_world_space.normalized()  # Ensure it's normalized
        print(f"Local Normal: {local_normal}, World Normal: {normal_in_world_space}")

        # Compute Theta (Rotation around Z)
        normal_xy = mathutils.Vector((normal_in_world_space.x, normal_in_world_space.y)).normalized()
        theta = math.atan2(normal_xy.y, normal_xy.x) + math.radians(90)
        theta = -self.normalize_angle(theta)
        print(f"Theta (Z-axis rotation): {math.degrees(theta)} degrees")

        # Apply theta rotation **to a copy of the normal vector**
        rotation_matrix_theta = self.rotation_matrix(theta, (0, 0, 1))
        rotated_normal = rotation_matrix_theta @ normal_in_world_space  # Rotate the normal only
        print(f"Rotated Normal (after theta): {rotated_normal}")

        # Compute Eta (Rotation around X)
        normal_yz = mathutils.Vector((rotated_normal.y, rotated_normal.z)).normalized()
        
        if normal_yz.length == 0:
            print("Warning: Normal is parallel to X-axis. Cannot calculate eta.")
            return None
        
        eta = math.atan2(normal_yz.x, normal_yz.y) + math.radians(-180)   
        eta = self.normalize_angle(eta)

        # Ensure Eta is between -90 and 0 degrees
        
        if bpy.context.scene.favor_negative_a and eta < math.radians(-90) or (eta > math.radians(45) and eta < math.radians(270)):
            print("Fixing positive eta...")
            theta = self.normalize_angle(theta + math.radians(180))  # Flip theta

            # Apply 180-degree rotation around Z to normal (not the object)
            rotation_matrix_theta = self.rotation_matrix(math.radians(180), (0, 0, 1))
            rotated_normal = rotation_matrix_theta @ rotated_normal  # Rotate the copied normal

            # Recalculate eta
            normal_yz = mathutils.Vector((rotated_normal.y, rotated_normal.z)).normalized()

            if normal_yz.length == 0:
                self.report({'ERROR'}, "Normal is parallel to X-axis after theta fix.")
                return None

            eta = math.atan2(normal_yz.x, normal_yz.y) + math.radians(-180)
            eta = self.normalize_angle(eta)
        elif not bpy.context.scene.favor_negative_a and eta > math.radians(90) or (eta < math.radians(-45) and eta > math.radians(-270)):
            print("Fixing negative eta...")
            theta = self.normalize_angle(theta - math.radians(180))  # Flip theta

            # Apply 180-degree rotation around Z to normal (not the object)
            rotation_matrix_theta = self.rotation_matrix(math.radians(180), (0, 0, 1))
            rotated_normal = rotation_matrix_theta @ rotated_normal  # Rotate the copied normal

            # Recalculate eta
            normal_yz = mathutils.Vector((rotated_normal.y, rotated_normal.z)).normalized()

            if normal_yz.length == 0:
                self.report({'ERROR'}, "Normal is parallel to X-axis after theta fix.")
                return None

            eta = math.atan2(normal_yz.x, normal_yz.y) + math.radians(-180)
            eta = self.normalize_angle(eta)
        # Apply final eta rotation **to normal vector only**
        rotation_matrix_eta = self.rotation_matrix(eta, (1, 0, 0))
        rotated_normal = rotation_matrix_eta @ rotated_normal  # Apply X-axis rotation

        # Final check on eta validity
        if math.radians(-90) <= eta <= math.radians(45):
            self.report({'INFO'}, f"{slicing_cube.name} rotations calulated")
            return theta, eta
        else:
            self.report({'ERROR'}, "A suitable rotation could not be found.")
            return 0, 0     
        
    def apply_boolean_intersect(self, original_mesh, slicing_cube):
        """
        Apply the Boolean Intersect operation on the given mesh using the slicing cube.
        Also applies the transformations (theta and eta) to the resulting sliced-off piece.
        """
        bpy.ops.object.select_all(action='DESELECT')
        original_mesh.select_set(True)
        bpy.context.view_layer.objects.active = original_mesh
        bool_intersect = original_mesh.modifiers.new(name="Intersect", type='BOOLEAN')
        bool_intersect.operation = 'INTERSECT'
        bool_intersect.object = slicing_cube
        bpy.context.view_layer.objects.active = original_mesh
        bpy.ops.object.modifier_apply(modifier=bool_intersect.name)
        
        self.report({'INFO'}, f"Boolean performed between {original_mesh} and {slicing_cube}")

        return None

    def apply_boolean_difference(self, selected_mesh, slicing_cube):
        """
        Apply the Boolean Difference operation on the selected mesh using the slicing cube.
        """
        bpy.context.view_layer.objects.active = selected_mesh
        bool_difference = selected_mesh.modifiers.new(name="Difference", type='BOOLEAN')
        bool_difference.operation = 'DIFFERENCE'
        bool_difference.object = slicing_cube
        bpy.ops.object.modifier_apply(modifier=bool_difference.name)

    def rotation_matrix(self, theta, axis):
        """
        Create a 4x4 rotation matrix that rotates by `theta` (radians) around `axis`.
        The axis must be a unit vector.
        """
        # Normalize the axis vector (make sure it’s a unit vector)
        axis = mathutils.Vector(axis).normalized()
        
        # Rodrigues' rotation formula components
        cos_theta = math.cos(theta)
        sin_theta = math.sin(theta)
        one_minus_cos = 1 - cos_theta
        
        # Elements of the rotation matrix
        ux, uy, uz = axis.x, axis.y, axis.z
        
        # Create the 4x4 rotation matrix
        rotation = mathutils.Matrix([
            [cos_theta + ux**2 * one_minus_cos, ux * uy * one_minus_cos - uz * sin_theta, ux * uz * one_minus_cos + uy * sin_theta, 0],
            [uy * ux * one_minus_cos + uz * sin_theta, cos_theta + uy**2 * one_minus_cos, uy * uz * one_minus_cos - ux * sin_theta, 0],
            [uz * ux * one_minus_cos - uy * sin_theta, uz * uy * one_minus_cos + ux * sin_theta, cos_theta + uz**2 * one_minus_cos, 0],
            [0, 0, 0, 1]
        ])
        
        return rotation
    
    def get_bbox_center(self, obj):
        if not obj:
            return None

        # Ensure Blender updates the object's data
        bpy.context.view_layer.update()

        # Convert bounding box (local space) to world space
        bbox_corners = [obj.matrix_world @ mathutils.Vector(corner) for corner in obj.bound_box]

        # Compute the center as a Blender Vector
        bbox_center = sum(bbox_corners, mathutils.Vector()) / 8  # Average of 8 points

        return bbox_center
    
    def get_set_bbox_center(self, obj):
        bbox_corners = obj.bound_box
        
        bbox_center = self.get_bbox_center(obj)
        
        bpy.context.scene.cursor.location = (bbox_center[0], bbox_center[1], bbox_center[2])
        self.report({'INFO'}, f"in order to set the bbox the cursos was moved to ({bbox_center[0]},{bbox_center[1]},{bbox_center[2]})")
        
        bpy.ops.object.select_all(action='DESELECT')
        obj.select_set(True)
        bpy.context.view_layer.objects.active = obj
        
        # Set the origin to the bounding box center
        bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
        z_coords = [obj.matrix_world @ v.co for v in obj.data.vertices]
        global_z_coords = [v.z for v in z_coords]
    
        # Find the minimum Z value in global coordinates
        min_z = min(global_z_coords)
    
        bpy.context.scene.cursor.location = (0, 0, bpy.context.scene.build_plate_distance)
       
        self.report({'INFO'}, f"{obj.name} set bbox to center returning {bbox_center[0]}, {bbox_center[1]}, {min_z}")
        
        return bbox_center[0], bbox_center[1], min_z
        
    def execute(self, context):
        scene = context.scene
        slicing_cubes = [bpy.data.objects[item.name] for item in scene.slicing_cubes_collection]
        
        selected_mesh = context.active_object
        selected_mesh.rotation_mode = 'YZX'
        self.set_origin_to_a_axis(selected_mesh)
        
        if not selected_mesh or selected_mesh.type != 'MESH':
            self.report({'ERROR'}, "Please select a valid mesh object!")
            return {'CANCELLED'}

        scene.selected_mesh = selected_mesh

        # Update cursor location directly
#        bpy.context.scene.cursor.location = (0, 0, scene.build_plate_distance)

        if not selected_mesh or selected_mesh.type != 'MESH':
            self.report({'ERROR'}, "Please select a mesh object to slice!")
            return {'CANCELLED'}
        
        if not slicing_cubes:
            self.report({'WARNING'}, "No slicing cubes found!")
            return {'CANCELLED'}
        
        if selected_mesh in slicing_cubes:
            self.report({'ERROR'}, "The selected object is a slicing cube. Please select a different object to slice.")
            return {'CANCELLED'}

        # Duplicate the selected mesh for the operation
        
        for slicing_cube in slicing_cubes:
            cube_data = next((item for item in scene.slicing_cubes_collection if item.name == slicing_cube.name), None)
            # Set the origin of the slicing cube to the 3D cursor
            
            if cube_data and cube_data.isProcessed:
                self.report({'INFO'}, f"Skipping already processed slicing cube: {slicing_cube.name}")
                continue  # Skip cubes that are already processed_lines
            
            mesh_copy = selected_mesh.copy()
            mesh_copy.data = selected_mesh.data.copy()
            scene.collection.objects.link(mesh_copy)
            
            self.report({'INFO'}, f"made copy: {mesh_copy.name}")
            
            slicing_cube.rotation_mode = 'YZX'
            self.set_origin_to_a_axis(slicing_cube)
         
            # Apply rotations based on theta and eta
            rotations = self.calculate_rotations(slicing_cube)
            if rotations is None:
                continue  # Skip if there is an issue with rotation calculation

            theta, eta = rotations
            
            
            
            self.set_origin_to_a_axis(mesh_copy)
            mesh_copy.rotation_mode = 'YZX'
            
            
            # Apply Boolean Intersect operation and apply transformations to the sliced-off piece
            self.apply_boolean_intersect(mesh_copy, slicing_cube)
            
            self.apply_rotations(mesh_copy, theta, eta)

            
            #apply rotations to sliced off part
            mesh_copy.location = (0,0,bpy.context.scene.build_plate_distance)
            x, y, min_z = self.get_set_bbox_center(mesh_copy)
        
            # Now, create a new sliced piece and assign the relevant information
            piece = scene.sliced_pieces_collection.add()
            piece.name = f"Sliced Piece {len(scene.sliced_pieces_collection)}"

            # Assign the sliced geometry (mesh_copy) to the piece
            piece.geometry_object = mesh_copy  # Store the geometry (mesh_copy) in the piece

            # Optionally assign metadata (like STL path)
            piece.stl_path = f"{mesh_copy.name}.stl"  # Placeholder for export path

            # Set the offsets based on the mesh's position
            piece.x_offset = x
            piece.y_offset = y
            piece.z_offset = min_z

            # Assign theta and eta (rotations) to the piece
            piece.theta = math.degrees(theta)  # Convert from radians to degrees for theta
            piece.eta = math.degrees(eta)      # Convert from radians to degrees for eta

            # Optionally, set default values for layer_height and infill_density
            piece.layer_height = 0.15 # Default layer height (modify as needed)
            piece.infill_density = 20.0  # Default infill density (modify as needed)

            # Apply Boolean Difference operation to the original mesh
            self.apply_boolean_difference(selected_mesh, slicing_cube)

            # Clear any transformations for the intersected mesh
            bpy.context.view_layer.objects.active = selected_mesh
            
            if cube_data:
                cube_data.isProcessed = True  # Flag this slicing cube as done
                self.report({'INFO'}, f"Slicing cube {slicing_cube.name} marked as processed.")

        # Cleanup by removing the duplicated original mesh
        
        #bounding box origin for selected mesh
        self.get_set_bbox_center(selected_mesh)
        
        self.report({'INFO'}, "Slicing completed successfully.")
        return {'FINISHED'}
    

###_______________________________________________________________
### GENERATE A
###________________________________________________________________

class SLICINGCUBE_OT_generate_gcode(bpy.types.Operator):
    bl_idname = "slicingcube.generate_gcode"
    bl_label = "Generate G-code"
    bl_description = "Generate G-code for the sliced model"
    bl_options = {'REGISTER', 'UNDO'}

    def get_bbox_center(self, obj):
        if not obj:
            return None

        # Ensure Blender updates the object's data
        bpy.context.view_layer.update()

        # Convert bounding box (local space) to world space
        bbox_corners = [obj.matrix_world @ mathutils.Vector(corner) for corner in obj.bound_box]

        # Compute the center as a Blender Vector
        bbox_center = sum(bbox_corners, mathutils.Vector()) / 8  # Average of 8 points

        return bbox_center
    
    def get_set_bbox_center(self, obj):
        bbox_corners = obj.bound_box
        
        bbox_center = self.get_bbox_center(obj)
        
        bpy.context.scene.cursor.location = (bbox_center[0], bbox_center[1], bbox_center[2])
        self.report({'INFO'}, f"in order to set the bbox the cursos was moved to ({bbox_center[0]},{bbox_center[1]},{bbox_center[2]})")
        
        bpy.ops.object.select_all(action='DESELECT')
        obj.select_set(True)
        bpy.context.view_layer.objects.active = obj
        
        # Set the origin to the bounding box center
        bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
        z_coords = [obj.matrix_world @ v.co for v in obj.data.vertices]
        global_z_coords = [v.z for v in z_coords]
    
        # Find the minimum Z value in global coordinates
        min_z = min(global_z_coords)
    
        bpy.context.scene.cursor.location = (0, 0, bpy.context.scene.build_plate_distance)
       
        self.report({'INFO'}, f"{obj.name} set bbox to center")
        
        return bbox_center[0], bbox_center[1], min_z
        

    def execute(self, context):
        scene = context.scene

        # Export the main model as an STL file
        stl_path = self.export_stl(scene.selected_mesh)
        if not stl_path:
            self.report({'ERROR'}, "Failed to export STL file for the main model.")
            return {'CANCELLED'}

        # Define the output G-code path for the main model
        output_gcode_path = os.path.join(os.path.dirname(stl_path), scene.selected_mesh.name + ".gcode")

        # Load the profile JSON
        profile_json = bpy.context.scene.fdmprinter_json_path

        # Generate G-code for the main model
        self.get_set_bbox_center(scene.selected_mesh)
        
        selected_mesh_location = scene.selected_mesh.location
        print(f"Selected mesh location: {selected_mesh_location}")

        self.slice_with_custom_coordinates(
            stl_path,
            output_gcode_path,
            selected_mesh_location.x,
            selected_mesh_location.y,
            0,  # Assuming Z-offset is zero for the main model
            profile_json,
            layer_height = 0.2,
            adhesion="brim",
            speed=30,
            support_enabled="true"
        )
        
        gcode_files = [output_gcode_path]
        
        self.report({'INFO'}, f"G-code generated for main model: {output_gcode_path}")

        # Iterate over sliced pieces and generate G-code for each
        for sliced_piece_item in list(reversed(scene.sliced_pieces_collection)):
            # Retrieve the Blender object associated with this slicing piece
            sliced_piece = sliced_piece_item.geometry_object
            
            if not sliced_piece:
                self.report({'WARNING'}, f"Could not find sliced piece object: {sliced_piece_item.name}")
                continue

            x, y, z = self.get_set_bbox_center(sliced_piece)
            
            #IMportant! ensures, no geometry is in teh -z region
            if sliced_piece_item.z_offset < 0:
                sliced_piece.location.z += -z
                
            # Export the sliced piece as an STL file
            stl_path_piece = self.export_stl(sliced_piece)
            if not stl_path_piece:
                self.report({'WARNING'}, f"Failed to export STL file for sliced piece: {sliced_piece_item.name}")
                continue

            # Define the output G-code path for the sliced piece
            output_gcode_piece = os.path.join(os.path.dirname(stl_path_piece), sliced_piece.name + ".gcode")

            # Get the sliced piece's location
            sliced_piece_location = sliced_piece.location
            print(f"Sliced piece location: {sliced_piece_location}")

            # Generate G-code for the sliced piece
            piece_gcode_path = self.slice_with_custom_coordinates(
                stl_path_piece,
                output_gcode_piece,
                x,
                y,
                sliced_piece_item.z_offset, 
                profile_json,
                layer_height = sliced_piece_item.layer_height,
                speed=sliced_piece_item.speed,
                adhesion="none",
                support_enabled="false"
            )
        
            # Append the safe turning position and rotation commands to the G-code file
            with open(piece_gcode_path, 'r') as file:
                gcode_lines = file.readlines()

            # Insert the safe turning position and rotation commands
            inital_comment = f"; THIS IS START OF {piece_gcode_path}\n"
            safe_position_command = "G1 F400 X0 Y0 Z140\n"
            rotation_command = f"G1 F400 A{sliced_piece_item.eta} {context.scene.rotary_axis_name}{sliced_piece_item.theta}\n"
            gcode_lines.insert(0, inital_comment)
            gcode_lines.insert(1, safe_position_command)
            gcode_lines.insert(2, rotation_command)

            # Write the modified G-code back to the file
            with open(output_gcode_piece, 'w') as file:
                file.writelines(gcode_lines)
            
            gcode_files.append(output_gcode_piece)
            

            self.report({'INFO'}, f"G-code generated for sliced piece: {output_gcode_piece}")
            
        # Define the path for the combined G-code file
        combined_gcode_path = os.path.join(os.path.dirname(stl_path), "combined_output.gcode")

        # Combine all G-code files into one (in proper sequential order)
        with open(combined_gcode_path, 'w') as combined_file:
            for tempfile in gcode_files:
                with open(tempfile,'r') as fi: combined_file.write(fi.read())
            # Add the "M2" command to signal the end of the program
            combined_file.write("M2\n")

        self.report({'INFO'}, f"Combined G-code file created: {combined_gcode_path}")
        
        return {'FINISHED'}

    def export_stl(self, obj):
        """
        Export the given object as an STL file.
        """
        print(obj)
        
        # Deselect all objects first
        bpy.ops.object.select_all(action='DESELECT')

        # Select the object to export
        obj.select_set(True)
        bpy.context.view_layer.objects.active = obj
        
        filename = obj.name + ".stl"

        # Define the export path
        stl_directory = os.path.dirname(bpy.context.scene.output_directory_path)
        stl_export_path = os.path.join(stl_directory, filename)

        # Export the STL file
        bpy.ops.wm.stl_export(filepath=stl_export_path, export_selected_objects=True)
        self.report({'INFO'}, f"STL file generated successfully: {stl_export_path}")    
    
        return stl_export_path
    
    

    def update_json_file(self, file_path, updates):
        """
        Update specific keys in a JSON file.

        :param file_path: Path to the JSON file to update.
        :param updates: Dictionary of key-value pairs to update in the JSON.
        """
        try:
            # Load the existing JSON data
            with open(file_path, 'r') as file:
                json_data = json.load(file)

            # Apply updates
            for key, value in updates.items():
                # Navigate to the nested key, if necessary
                keys = key.split('.')
                data = json_data
                for sub_key in keys[:-1]:
                    data = data.get(sub_key, {})
                # Update the final key
                data[keys[-1]] = value

            # Save the updated JSON back to the file
            with open(file_path, 'w') as file:
                json.dump(json_data, file, indent=4)

            self.report({'INFO'}, f"Successfully updated {file_path}. Updates: {updates}")
        except Exception as e:
            print(f"Error updating JSON file: {e}")

    
    def slice_with_custom_coordinates(self, input_stl, output_gcode, x_offset, y_offset, z_offset, profile_json, layer_height=0.225, speed=40, infill=15, adhesion="brim", support_enabled=False):
        """
        Slice the STL file and generate G-code using CuraEngine.
        """
        
        # Dictionary of updates you want to make
        updates = {
            "settings.command_line_settings.children.mesh_position_x.default_value": 0,
            "settings.command_line_settings.children.mesh_position_y.default_value": 0,
            "settings.command_line_settings.children.mesh_position_z.default_value": 0,
            "settings.speed.children.speed_print.children.default_value": f"{speed:2f}",
            "settings.platform_adhesion.children.adhesion_type.default_value": adhesion,
            "settings.resolution.children.layer_height.default_value": f"{layer_height:.2f}",
            "settings.support.children.support_enable.default_value": support_enabled
            
        }
        
        temp_path = input_stl + ".gcode"

        self.update_json_file(bpy.context.scene.fdmprinter_json_path, updates)
        
        command = [
            bpy.context.scene.cura_executable_path,
            "slice",
            "-p",
            "-s", 'machine_start_gcode=',   # This line disables start g-code
            "-s", 'machine_end_gcode=M2',
            "-s", 'roofing_layer_count=3',
            "-s", 'adhesion_type=none',
#            "-s", 'brim_width=0',
#            "-s", 'layer_height_0=0.2',
#            "-s", 'layer_height=0.2',
#            "-s", f'speed_print={speed}',
#            "-s", 'speed_travel=60',
            "-s", 'machine_extruder_count=1',
            "-s", 'mesh_rotation_matrix="[[1,0,0],[0,1,0],[0,0,1]]"',
            "-s", 'machine_center_is_zero=true',
            "-s", 'mesh_position_x=0',
            "-s", 'mesh_position_y=0',
            "-s", 'mesh_position_z=0',
            "-s", 'center_object=true',
            "-s", 'material_shrinkage_percentage_xy=0',
            "-s", 'material_shrinkage_percentage_z=0',
            "-s", 'infill_mesh=false',
            "-s", 'cutting_mesh=false',
            "-s", 'anti_overhang_mesh=false',
            "-s", 'machine_width=220',
            "-s", 'machine_depth=220',
            "-s", 'machine_height=220',
            "-s", 'machine_heated_bed=true',
            "-s", 'machine_heated_build_volume=false',
            "-s", 'material_guid=false',
            "-s", 'acceleration_enabled=false',
            "-s", 'jerk_enabled=false',
            "-s", 'relative_extrusion=false',
            "-s", 'machine_scale_fan_speed_zero_to_one=1',
            "-s", 'machine_nozzle_temp_enabled=false',
            "-s", 'adaptive_layer_height_enabled=false',
            "-s", 'slicing_tolerance=0.05',
            "-s", 'support_enable=false',
            "-s", 'support_mesh=false',
            "-s", 'remove_empty_first_layers=false',
            "-s", 'skirt_brim_extruder_nr=1',
            "-s", 'prime_tower_brim_enable=false',
            "-s", 'machine_gcode_flavor=Mach3',
            "-s", 'machine_use_extruder_offset_to_offset_coords=false',
            "-s", 'ppr_enable=false',
            "-s", 'retraction_prime_speed=10',
            "-s", 'machine_firmware_retract=false',
            "-s", 'machine_name=JAiro',
            "-s", 'machine_nozzle_offset_x=0',
            "-s", 'machine_nozzle_offset_y=0',
            "-s", 'machine_nozzle_offset_z=0',
            "-s", 'machine_always_write_active_tool=true',
            "-s", 'machine_max_feedrate_x=3600',
            "-s", 'machine_max_feedrate_y=3600',
            "-s", 'machine_max_feedrate_z=3600',
            "-s", 'machine_max_feedrate_e=3600',
            "-s", 'machine_max_acceleration_x=500',
            "-s", 'machine_max_acceleration_y=500',
            "-s", 'machine_max_acceleration_z=500',
            "-s", 'machine_max_acceleration_e=500',
            "-s", 'machine_max_jerk_xy=20',
            "-s", 'machine_max_jerk_z=20',
            "-s", 'machine_max_jerk_e=20',
            "-s", 'machine_minimum_feedrate=0',
            "-s", 'machine_acceleration=1000',
            "-s", 'machine_start_gcode_first=""',
            "-s", 'material_bed_temp_prepend=100',
            "-s", 'material_print_temperature_layer_0=200',
            "-s", 'material_bed_temperature_layer_0=60',
            "-s", 'cool_min_layer_time=10',
            "-s", 'cool_min_layer_time_fan_speed_max=100',
            "-s", 'cool_fan_speed_0=100',
            "-s", 'cool_fan_speed_min=0',
            "-s", 'cool_fan_speed_max=0',
            "-s", 'cool_min_speed=10',
            "-s", 'cool_fan_full_layer=2',
            "-s", 'cool_fan_enabled=true',
            "-s", 'retraction_enable=true',
            "-s", 'retraction_amount=3',
            "-s", 'retraction_extra_prime_amount=0',
            "-s", 'retraction_retract_speed=5',
            "-s", 'retraction_hop=true',
            "-s", 'retraction_min_travel=3',
            "-s", 'retraction_extrusion_window=4.5',
            "-s", 'retraction_count_max=25',
            "-s", 'retraction_hop_after_extruder_switch=true',
            "-s", 'switch_extruder_extra_prime_amount=0',
            "-s", 'switch_extruder_retraction_amount=0',
            "-s", 'switch_extruder_retraction_speed=5',
            "-s", 'switch_extruder_prime_speed=50',
            "-s", 'retraction_hop_after_extruder_switch_height=0',
            "-s", 'wipe_retraction_enable=false',
            "-s", 'wipe_retraction_amount=5',
            "-s", 'wipe_retraction_retract_speed=10',
            "-s", 'wipe_retraction_prime_speed=10',
            "-s", 'wipe_retraction_extra_prime_amount=5',
            "-s", 'wipe_pause=0',
            "-s", 'wipe_hop_enable=false',
            "-s", 'wipe_hop_amount=0',
            "-s", 'wipe_hop_speed=0',
            "-s", 'wipe_brush_pos_x=0',
            "-s", 'wipe_brush_pos_y=0',
            "-s", 'wipe_brush_pos_z=0',
            "-s", 'wipe_repeat_count=0',
            "-s", 'wipe_move_distance=0',
            "-s", 'max_extrusion_before_wipe=10',
            "-s", 'clean_between_layers=false',
            "-s", 'material_bed_temp_wait=40',
            "-s", 'material_print_temp_prepend=200',
            "-s", 'material_print_temp_wait=200',
            "-s", 'machine_extruders_share_nozzle=false',
            "-s", 'support_infill_extruder_nr=0',
            "-s", 'support_infill_angles="[]"',
            "-s", 'support_extruder_nr_layer_0=0',
            "-s", 'support_roof_extruder_nr=0',
            "-s", 'support_roof_angles="[]"',
            "-s", 'support_roof_pattern=lines',
            "-s", 'support_bottom_extruder_nr=0',
            "-s", 'support_bottom_pattern=lines',
            "-s", 'support_bottom_angles="[]"',
            "-s", 'layer_start_x=0',
            "-s", 'layer_start_y=0',
            "-s", 'layer_start_z=0',
            "-s", 'prime_blob_enable=false',
            "-s", 'print_sequence=all_at_once',
            "-s", 'raft_surface_extruder_nr=0',
            "-s", 'magic_spiralize=false',
            "-s", 'extruder_nr=0',
            "-s", 'material_diameter=1.75',
            "-s", 'machine_extruder_cooling_fan_number=0',
            "-s", 'machine_extruder_prestart_code=',
            "-s", 'machine_extruder_start_code=',
            "-s", 'machine_extruder_end_code=',
            "-s", 'machine_extruder_start_code_duration=0',
            "-s", 'machine_extruder_change_duration=0',
            "-s", 'machine_extruder_start_pos_x=0',
            "-s", 'machine_extruder_start_pos_y=0',
            "-s", 'machine_extruder_start_pos_z=0',
            "-s", 'reset_flow_duration=0',
            "-s", 'support_z_seam_away_from_model=true',
            "-s", 'min_wall_line_width="0.3"',
            "-j", profile_json,
            "-l", input_stl,
            "-o", temp_path
        ]
        
        print("Generated CuraEngine Command:", " ".join(command))
        
        # Run the command
        result = subprocess.run(command, capture_output=True, text=True)

        # Print the output and error (if any)
        # Print the output and error (if any)
        print("Command:", " ".join(command))
        print("Output:", result.stdout)
        print("Error:", result.stderr)

        if result.returncode != 0:
            self.report({'ERROR'}, "Error while slicing the model. Check the console for details.")
            return None

        # Remove the first few lines of the G-code (if needed)
        self.clean_gcode(temp_path, output_gcode, z_offset=z_offset, extruder_to_c=True, lines_to_remove=22, lines_to_remove_end=6)
        
        return output_gcode
    
    def remove_unwanted_commands(self, line):
        """
        Remove unwanted commands (e.g., M107, M106, M104) from the G-code line.
        """
        if line.startswith("M140") or line.startswith("M82") or line.startswith("M107") or line.startswith("M106") or line.startswith("M104"):
            return True
        return False


    def apply_z_offset(self, line, z_offset):
        """
        Apply the Z-offset to all Z-axis movements in the G-code line.
        """
        if " Z" in line:
            parts = line.split()
            for i, part in enumerate(parts):
                if part.startswith("Z"):
                    try:
                        z_value = float(part[1:])
                        z_value += z_offset  # Apply Z offset
                        parts[i] = f"Z{z_value:.3f}"
                    except ValueError:
                        pass
            line = " ".join(parts) + "\n"
        return line


    def replace_a_with_c(self, line):
        """
        Replace the A-axis with C-axis in the G-code line.
        """
        if "E" in line:  # Check for "A" in the line (without worrying about spaces)
            parts = line.split()
            for i, part in enumerate(parts):
                if part.startswith("E"):  # Check for "A" at the start of the part
                    try:
                        a_value = float(part[1:])  # Extract the numeric value after "A"
                        parts[i] = f"{bpy.context.scene.extruder_axis_name}{a_value:.3f}"  # Replace "A" with "C"
                    except ValueError:
                        pass
            line = " ".join(parts) + "\n"  # Rebuild the line
        return line
    
    def clean_gcode(self, input_gcode, output_gcode, z_offset, extruder_to_c=True, lines_to_remove=0, lines_to_remove_end=0):
        """
        Main function to clean up G-code by applying transformations such as Z-offset, replacing axes, and removing unwanted commands.
        
        :param input_gcode: Path to the input G-code file.
        :param output_gcode: Path to the output G-code file.
        :param z_offset: Float value to offset Z-axis positions.
        :param extruder_to_c: Boolean to replace A-axis with C-axis.
        :param lines_to_remove: Number of lines to remove from the end of the G-code.
        """
        # First, read all lines from the input G-code into memory
        with open(input_gcode, 'r') as infile:
            lines = infile.readlines()

        # Process all lines sequentially
        processed_lines = []
        for line in lines:
            # Step 1: Remove unwanted commands
            if self.remove_unwanted_commands(line):
                continue  # Skip lines that should be removed

            # Step 2: Apply Z-offset to Z-movements
            line = self.apply_z_offset(line, z_offset)

            # Step 3: Replace A-axis with C-axis (extruder -> rotary axis)
            if extruder_to_c:
                line = self.replace_a_with_c(line)

            # Add the modified line to the processed list
            processed_lines.append(line)
        
        self.report({'INFO'}, f"Adjusted Z-offser by {z_offset}")
        
        # Step 4: Optionally, remove the first 'lines_to_remove' lines from the processed lines
        if lines_to_remove > 0 and len(processed_lines) > lines_to_remove:
            processed_lines = processed_lines[lines_to_remove:]

        # Step 5: Optionally, remove the last 'lines_to_remove' lines from the processed lines
        if lines_to_remove_end > 0 and len(processed_lines) > lines_to_remove_end:
            processed_lines = processed_lines[:-lines_to_remove_end]

        # Write the processed lines back to the output file
        with open(output_gcode, 'w') as outfile:
            outfile.writelines(processed_lines)

        self.report({'INFO'}, f"Done adjusting G-code file: {output_gcode}")
# -------------------------------------------------------------------
# (6) sliced off parts
# -------------------------------------------------------------------
# Define a data structure for sliced pieces
class SlicedPieceItem(bpy.types.PropertyGroup):
    name: bpy.props.StringProperty(name="Name", default="Sliced Piece")
    stl_path: bpy.props.StringProperty(name="STL Path", default="")
    theta: bpy.props.FloatProperty(name="Theta (Z Rotation)", default=0.0)
    eta: bpy.props.FloatProperty(name="Eta (X Rotation)", default=0.0)
    x_offset: bpy.props.FloatProperty(name="X Offset", default=0.0)
    y_offset: bpy.props.FloatProperty(name="Y Offset", default=0.0)
    z_offset: bpy.props.FloatProperty(name="Z Offset", default=0.0)
    geometry_object: bpy.props.PointerProperty(
        name="Geometry Object",
        type=bpy.types.Object,
        description="Reference to the sliced geometry object"
    )

    # Add the missing properties
    layer_height: bpy.props.FloatProperty(
        name="Layer Height",
        description="The layer height for 3D printing",
        default=0.2,  # Default value in mm
        min=0.05,  # Minimum layer height (example limit)
        max=1.0    # Maximum layer height (example limit)
    )
    speed: bpy.props.FloatProperty(
        name="Speed",
        description="Print speed",
        default=60,  # Default value in mm
        min=0,  # Minimum layer height (example limit)
        max=1100    # Maximum layer height (example limit)
    )
    infill_density: bpy.props.FloatProperty(
        name="Infill Percentage",
        description="Infill density as a percentage",
        default=20.0,  # Default value
        min=0.0,       # 0% infill
        max=100.0      # 100% solid
    )
    
class SLICEDPIECE_UL_items(bpy.types.UIList):
    def draw_item(self, context, layout, data, item, icon, active_data, active_propname, index):
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            layout.label(text=item.name, icon='MESH_CUBE')
        elif self.layout_type in {'GRID'}:
            layout.alignment = 'CENTER'
            layout.label(text="", icon='MESH_CUBE')

class SLICEDPIECE_OT_add(bpy.types.Operator):
    bl_idname = "slicedpiece.add"
    bl_label = "Add Sliced Piece"
    bl_description = "Add a new sliced piece to the collection"

    def execute(self, context):
        scene = context.scene
        # Get the selected object (the mesh you want to slice)
        selected_mesh = context.active_object
        
        # Ensure the selected object is a mesh
        if not selected_mesh or selected_mesh.type != 'MESH':
            self.report({'ERROR'}, "Please select a mesh object to slice!")
            return {'CANCELLED'}

        # Create a new SlicedPieceItem and add it to the collection
        piece = scene.sliced_pieces_collection.add()
        scene.sliced_pieces_index = len(scene.sliced_pieces_collection) - 1

        return {'FINISHED'}
    
class SLICEDPIECE_OT_remove(bpy.types.Operator):
    bl_idname = "slicedpiece.remove"
    bl_label = "Remove Sliced Piece"
    bl_description = "Remove the selected sliced piece from the collection"

    def execute(self, context):
        scene = context.scene
        if scene.sliced_pieces_index >= 0:
            scene.sliced_pieces_collection.remove(scene.sliced_pieces_index)
            scene.sliced_pieces_index = min(scene.sliced_pieces_index, len(scene.sliced_pieces_collection) - 1)
        return {'FINISHED'}
    
# -------------------------------------------------------------------
# (6) Main Panel
# -------------------------------------------------------------------
class SLICINGCUBE_OT_reset_slicing(bpy.types.Operator):
    """Reset slicing cubes and sliced parts"""
    bl_idname = "slicingcube.reset_slicing"
    bl_label = "Reset Slicing Data"

    def execute(self, context):
        self.reset_slicing_data(context)  # Call the reset function here
        return {'FINISHED'}
    
    def reset_slicing_data(self, context):
        """
        Resets the slicing cubes and sliced pieces in the scene.
        """
        scene = context.scene
        
        # Clear the slicing cubes collection
        slicing_cubes_collection = scene.slicing_cubes_collection
        for slicing_cube_item in slicing_cubes_collection:
            # Get the actual Blender object by name
            slicing_cube_object = bpy.data.objects.get(slicing_cube_item.name)
            if slicing_cube_object and slicing_cube_object.users_collection:  # Check if the object exists and is in any collection
                # Remove from the scene collection
                for coll in slicing_cube_object.users_collection:
                    coll.objects.unlink(slicing_cube_object)
        
        # Clear the sliced pieces collection
        sliced_pieces_collection = scene.sliced_pieces_collection
        for sliced_piece_item in sliced_pieces_collection:
            # Get the actual Blender object by name
            sliced_piece_object = bpy.data.objects.get(sliced_piece_item.name)
            if sliced_piece_object and sliced_piece_object.users_collection:  # Check if the object exists and is in any collection
                # Remove from the scene collection
                for coll in sliced_piece_object.users_collection:
                    coll.objects.unlink(sliced_piece_object)

        # Clear the custom collections themselves
        scene.slicing_cubes_collection.clear()
        scene.sliced_pieces_collection.clear()

        self.report({'INFO'}, "Slicing data reset successfully.")
        return {'FINISHED'}
    
class VIEW3D_PT_5AxisPrinterSetup(bpy.types.Panel):
    """Main panel for 5-Axis Printer, Slicing Cubes, and Slice"""
    bl_label = "5-Axis Printer INdexed Slicer"
    bl_idname = "VIEW3D_PT_5axis_3dprinter_bottom_face_down"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "5-Axis Printer"

    def draw(self, context):
        layout = self.layout
        scene = context.scene
        
        layout.prop(scene, "rotary_axis_name", text="2nd Rot Letter")
        layout.prop(scene, "extruder_axis_name", text="Extruder Letter")
        layout.prop(scene, "favor_negative_a")
        
        # Add the distance parameter as an input field
        layout.label(text="Set the Distance to Build Plate (A-axis):")
        layout.prop(scene, "build_plate_distance", text="Distance")

        # Section: Create 5-Axis Setup
        layout.label(text="Create 5-Axis Setup")
        
        # Section: Slicing Cubes
        layout.separator()
        layout.label(text="Slicing Cubes")

        row = layout.row()
        row.template_list(
            "SLICINGCUBE_UL_items", "",
            scene, "slicing_cubes_collection",
            scene, "slicing_cubes_index",
            rows=3
        )
        col = row.column(align=True)
        col.operator("slicingcube.add", icon='ADD', text="")
        col.operator("slicingcube.remove", icon='REMOVE', text="")

        # Section: Calculate Theta
        layout.separator()
        layout.label(text="Calculate Rotations")
        layout.operator("slicingcube.calculate_theta", text="Calculate Rotations", icon="DRIVER_ROTATIONAL_DIFFERENCE")

        # Section: Sliced Pieces
        layout.separator()
        layout.label(text="Sliced Pieces")

        # Display sliced pieces
        row = layout.row()
        row.template_list(
            "SLICEDPIECE_UL_items", "",
            scene, "sliced_pieces_collection",
            scene, "sliced_pieces_index",
            rows=3
        )
        col = row.column(align=True)
        col.operator("slicedpiece.add", icon='ADD', text="")
        col.operator("slicedpiece.remove", icon='REMOVE', text="")

        # Display properties of the selected sliced piece
        if scene.sliced_pieces_index >= 0 and len(scene.sliced_pieces_collection) > 0:
            piece = scene.sliced_pieces_collection[scene.sliced_pieces_index]
           
                # Display properties of the selected sliced piece
            layout.prop(piece, "geometry_object", text="Geometry Object")
            layout.prop(piece, "layer_height")
            layout.prop(piece, "speed")
            layout.prop(piece, "infill_density")
            layout.prop(piece, "theta")
            layout.prop(piece, "eta")
            layout.prop(piece, "x_offset")
            layout.prop(piece, "y_offset")
            layout.prop(piece, "z_offset")

        # Section: Generate G-code
        layout.separator()
        layout.prop(scene, "cura_executable_path")
        layout.prop(scene, "fdmprinter_json_path")
        layout.prop(scene, "output_directory_path")
        layout.label(text="Generate G-code")
#        layout.prop(
        layout.operator("slicingcube.generate_gcode", text="Generate G-code", icon="FILE_SCRIPT")
        
        # Reset Button
        layout.separator()
        layout.label(text="Reset Slicing Data")
        layout.operator("slicingcube.reset_slicing", text="Reset Slicing Data", icon="PANEL_CLOSE")

# -------------------------------------------------------------------
# Registration
# -------------------------------------------------------------------
classes = (
    SlicingCubeItem,
    SLICINGCUBE_UL_items,
    SLICINGCUBE_OT_add,
    SLICINGCUBE_OT_remove,
    SLICINGCUBE_OT_slice,
    SLICINGCUBE_OT_generate_gcode,
    VIEW3D_PT_5AxisPrinterSetup,
    SlicedPieceItem,
    SLICEDPIECE_OT_add,
    SLICEDPIECE_OT_remove,
    SLICEDPIECE_UL_items,
    SLICINGCUBE_OT_reset_slicing
    
)

def register():
    for cls in classes:
        bpy.utils.register_class(cls)

    bpy.types.Scene.slicing_cubes_collection = bpy.props.CollectionProperty(type=SlicingCubeItem)
    bpy.types.Scene.slicing_cubes_index = bpy.props.IntProperty(default=-1)
    bpy.types.Scene.sliced_pieces_collection = CollectionProperty(type=SlicedPieceItem)
    bpy.types.Scene.sliced_pieces_index = bpy.props.IntProperty(default=0)
    
    bpy.types.Scene.build_plate_distance = bpy.props.FloatProperty(
        name="Build Plate Distance",
        default=-35.8,  # Default distance
        min=-250,  # Minimum value (can adjust if needed)
        max=250,  # Maximum value (can adjust as per needs)
        unit='LENGTH',  # Unit is length in this case (millimeters)
        description="Distance between the build plate and the A-axis"
    )
    bpy.types.Scene.selected_mesh = bpy.props.PointerProperty(
        name="Selected Mesh",
        type=bpy.types.Object,
        description="The mesh selected for slicing"
    )
    # Register the scene properties for global parameters
    bpy.types.Scene.x_width = bpy.props.FloatProperty(
        name="X Width",
        description="Width of the build plate in mm",
        default=220.0,  # Set a default value (example: 200mm)
        min=0.1,        # Optional: Add a minimum limit
        max=1000.0      # Optional: Add a maximum limit
    )

    bpy.types.Scene.y_depth = bpy.props.FloatProperty(
        name="Y Depth",
        description="Depth of the build plate in mm",
        default=220.0,  # Set a default value (example: 200mm)
        min=0.1,        # Optional: Add a minimum limit
        max=1000.0      # Optional: Add a maximum limit
    )
    bpy.types.Scene.rotary_axis_name = bpy.props.StringProperty(
        name="B Axis Name",
        description="The name of the axis to replace the extruder (default: C)",
        default="B"
    )
    bpy.types.Scene.extruder_axis_name = bpy.props.StringProperty(
        name="Extruder Axis Name",
        description="The name of the axis to replace the extruder (default: C)",
        default="E"
    )
    bpy.types.Scene.favor_negative_a = bpy.props.BoolProperty(
        name="Favor negative A rotations",
        description="Which way to turn A",
        default = True,
        update=lambda self, context: print(f"Updated: {self.favor_negative_a}")
    )
    bpy.types.Scene.cura_executable_path = bpy.props.StringProperty(
        name="Cura executable file path",
        description="Choose a file",
        default="/Users/jairo/CuraEngine/build/Release/CuraEngine",
        subtype='FILE_PATH'
    )
    bpy.types.Scene.fdmprinter_json_path = bpy.props.StringProperty(
        name="fdmprinter JSON file path",
        description="Choose a file",
        default="/Users/jairo/Documents/4_5th_axis/Addon/fdmprinter.def.json",
        subtype='FILE_PATH'
    )
    bpy.types.Scene.output_directory_path = bpy.props.StringProperty(
        name="Output directory path",
        description="Choose a file",
        default="",
        subtype='DIR_PATH'
    )

    
def unregister():
    # Unregister collection properties
    del bpy.types.Scene.slicing_cubes_collection
    del bpy.types.Scene.slicing_cubes_index
    del bpy.types.Scene.sliced_pieces_collection
    del bpy.types.Scene.sliced_pieces_index

    # Unregister custom properties
    del bpy.types.Scene.build_plate_distance
    del bpy.types.Scene.selected_mesh
    del bpy.types.Scene.x_width
    del bpy.types.Scene.y_depth
    del bpy.types.Scene.c_axis_name
    del bpy.types.Scene.extruder_axis_name
    del bpy.types.Scene.favor_negative_a
    
    # Unregister classes
    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)

if __name__ == "__main__":
    register()
    
 
    