#include "vr_cobotics.h"

#include <cgv/signal/rebind.h>
#include <cgv/base/register.h>
#include <cgv/math/ftransform.h>
#include <cgv/utils/scan.h>
#include <cgv/utils/options.h>
#include <cgv/utils/file.h>
#include <cgv/gui/dialog.h>
#include <cgv/gui/file_dialog.h>
#include <cgv/render/attribute_array_binding.h>
#include <cgv_gl/sphere_renderer.h>
#include <cgv/media/mesh/simple_mesh.h>
#include <cg_vr/vr_events.h>

#include <random>
#include <sstream>

#include "intersection.h"

template <typename T>
double time_stamp(const T& t_start) {
	auto now = std::chrono::steady_clock::now();
	return static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(now - t_start).count()) / 1e6;
}

bool hasEnding(std::string const &str, std::string const &ending) {
	if (str.length() >= ending.length()) {
		return str.compare(str.length() - ending.length(), ending.length(), ending) == 0;
	}
	return false;
}

void vr_cobotics::change_box_extents(Axis axis,int ci) {
	for (size_t i = 0; i < intersection_points.size(); ++i) {
		if (intersection_controller_indices[i] != ci)
			continue;
		// extract box index
		int bi = intersection_box_indices[i];
		box3 b = movable_boxes[bi];
		float extent = b.get_max_pnt()[axis] - b.get_min_pnt()[axis];
		float center = b.get_center()[axis];

		const float size_limit = 0.5;
		float new_extent = std::max(std::fmod(extent + edit_box_step, size_limit), edit_box_step);
		movable_boxes[bi].ref_max_pnt()[axis] = center + new_extent * 0.5f;
		movable_boxes[bi].ref_min_pnt()[axis] = center - new_extent * 0.5f;
		// update intersection points
		//intersection_points[i] = rotation * (intersection_points[i] - last_pos) + pos;
		break; //change only the first box
	}
}

void vr_cobotics::delete_box(int bi)
{
	movable_boxes.erase(movable_boxes.begin() + bi);
	movable_box_colors.erase(movable_box_colors.begin() + bi);
	movable_box_translations.erase(movable_box_translations.begin() + bi);
	movable_box_rotations.erase(movable_box_rotations.begin() + bi);
}

size_t vr_cobotics::clear_intersections(int ci)
{
	size_t i = 0;
	while (i < intersection_points.size()) {
		if (intersection_controller_indices[i] == ci) {
			intersection_points.erase(intersection_points.begin() + i);
			intersection_colors.erase(intersection_colors.begin() + i);
			intersection_box_indices.erase(intersection_box_indices.begin() + i);
			intersection_controller_indices.erase(intersection_controller_indices.begin() + i);
		}
		else
			++i;
	}
	return i;
}

void vr_cobotics::init_cameras(vr::vr_kit* kit_ptr)
{
	vr::vr_camera* camera_ptr = kit_ptr->get_camera();
	if (!camera_ptr)
		return;
	nr_cameras = camera_ptr->get_nr_cameras();
	frame_split = camera_ptr->get_frame_split();
	for (int i = 0; i < nr_cameras; ++i) {
		std::cout << "camera " << i << "(" << nr_cameras << "):" << std::endl;
		camera_ptr->put_camera_intrinsics(i, false, &focal_lengths[i](0), &camera_centers[i](0));
		camera_ptr->put_camera_intrinsics(i, true, &focal_lengths[2 + i](0), &camera_centers[2 + i](0));
		std::cout << "  fx=" << focal_lengths[i][0] << ", fy=" << focal_lengths[i][1] << ", center=[" << camera_centers[i] << "]" << std::endl;
		std::cout << "  fx=" << focal_lengths[2+i][0] << ", fy=" << focal_lengths[2+i][1] << ", center=[" << camera_centers[2+i] << "]" << std::endl;
		float camera_to_head[12];
		camera_ptr->put_camera_to_head_matrix(i, camera_to_head);
		kit_ptr->put_eye_to_head_matrix(i, camera_to_head);
		camera_to_head_matrix[i] = vr::get_mat4_from_pose(camera_to_head);
		std::cout << "  C2H=" << camera_to_head_matrix[i] << std::endl;
		camera_ptr->put_projection_matrix(i, false, 0.001f, 10.0f, &camera_projection_matrix[i](0, 0));
		camera_ptr->put_projection_matrix(i, true, 0.001f, 10.0f, &camera_projection_matrix[2+i](0, 0));
		std::cout << "  dP=" << camera_projection_matrix[i] << std::endl;
		std::cout << "  uP=" << camera_projection_matrix[2+i] << std::endl;
	}
	post_recreate_gui();
}

void vr_cobotics::start_camera()
{
	if (!vr_view_ptr)
		return;
	vr::vr_kit* kit_ptr = vr_view_ptr->get_current_vr_kit();
	if (!kit_ptr)
		return;
	vr::vr_camera* camera_ptr = kit_ptr->get_camera();
	if (!camera_ptr)
		return;
	if (!camera_ptr->start())
		cgv::gui::message(camera_ptr->get_last_error());
}

void vr_cobotics::stop_camera()
{
	if (!vr_view_ptr)
		return;
	vr::vr_kit* kit_ptr = vr_view_ptr->get_current_vr_kit();
	if (!kit_ptr)
		return;
	vr::vr_camera* camera_ptr = kit_ptr->get_camera();
	if (!camera_ptr)
		return;
	if (!camera_ptr->stop())
		cgv::gui::message(camera_ptr->get_last_error());
}

/// compute intersection points of controller ray with movable boxes
void vr_cobotics::compute_intersections(const vec3& origin, const vec3& direction, int ci, const rgb& color)
{
	for (size_t i = 0; i < movable_boxes.size(); ++i) {
		vec3 origin_box_i = origin - movable_box_translations[i];
		movable_box_rotations[i].inverse_rotate(origin_box_i);
		vec3 direction_box_i = direction;
		movable_box_rotations[i].inverse_rotate(direction_box_i);
		float t_result;
		vec3  p_result;
		vec3  n_result;
		if (cgv::media::ray_axis_aligned_box_intersection(
			origin_box_i, direction_box_i,
			movable_boxes[i],
			t_result, p_result, n_result, 0.000001f)) {

			// transform result back to world coordinates
			movable_box_rotations[i].rotate(p_result);
			p_result += movable_box_translations[i];
			movable_box_rotations[i].rotate(n_result);

			// store intersection information
			intersection_points.push_back(p_result);
			intersection_colors.push_back(color);
			intersection_box_indices.push_back((int)i);
			intersection_controller_indices.push_back(ci);
		}
	}
}

/// register on device change events
void vr_cobotics::on_device_change(void* kit_handle, bool attach)
{
	if (attach) {
		if (last_kit_handle == 0) {
			vr::vr_kit* kit_ptr = vr::get_vr_kit(kit_handle);
			init_cameras(kit_ptr);
			if (kit_ptr) {
				last_kit_handle = kit_handle;
				left_deadzone_and_precision = kit_ptr->get_controller_throttles_and_sticks_deadzone_and_precision(0);
				cgv::gui::ref_vr_server().provide_controller_throttles_and_sticks_deadzone_and_precision(kit_handle, 0, &left_deadzone_and_precision);
				post_recreate_gui();
			}
		}
	}
	else {
		if (kit_handle == last_kit_handle) {
			last_kit_handle = 0;
			post_recreate_gui();
		}
	}
}

/// construct boxes that represent a table of dimensions tw,td,th and leg width tW
void vr_cobotics::construct_table(float tw, float td, float th, float tW) {
	// construct table
	rgb table_clr(0.3f, 0.2f, 0.0f);
	boxes.push_back(box3(
		vec3(-0.5f*tw - 2 * tW, th - tW, -0.5f*td - 2 * tW),
		vec3(0.5f*tw + 2 * tW, th, 0.5f*td + 2 * tW)));
	box_colors.push_back(table_clr);

	boxes.push_back(box3(vec3(-0.5f*tw, 0, -0.5f*td), vec3(-0.5f*tw - tW, th - tW, -0.5f*td - tW)));
	boxes.push_back(box3(vec3(-0.5f*tw, 0, 0.5f*td), vec3(-0.5f*tw - tW, th - tW, 0.5f*td + tW)));
	boxes.push_back(box3(vec3(0.5f*tw, 0, -0.5f*td), vec3(0.5f*tw + tW, th - tW, -0.5f*td - tW)));
	boxes.push_back(box3(vec3(0.5f*tw, 0, 0.5f*td), vec3(0.5f*tw + tW, th - tW, 0.5f*td + tW)));
	box_colors.push_back(table_clr);
	box_colors.push_back(table_clr);
	box_colors.push_back(table_clr);
	box_colors.push_back(table_clr);
}

/// construct boxes that represent a room of dimensions w,d,h and wall width W
void vr_cobotics::construct_room(float w, float d, float h, float W, bool walls, bool ceiling) {
	// construct floor
	boxes.push_back(box3(vec3(-0.5f*w, -W, -0.5f*d), vec3(0.5f*w, 0, 0.5f*d)));
	box_colors.push_back(rgb(0.2f, 0.2f, 0.2f));

	if(walls) {
		// construct walls
		boxes.push_back(box3(vec3(-0.5f*w, -W, -0.5f*d - W), vec3(0.5f*w, h, -0.5f*d)));
		box_colors.push_back(rgb(0.8f, 0.5f, 0.5f));
		boxes.push_back(box3(vec3(-0.5f*w, -W, 0.5f*d), vec3(0.5f*w, h, 0.5f*d + W)));
		box_colors.push_back(rgb(0.8f, 0.5f, 0.5f));

		boxes.push_back(box3(vec3(0.5f*w, -W, -0.5f*d - W), vec3(0.5f*w + W, h, 0.5f*d + W)));
		box_colors.push_back(rgb(0.5f, 0.8f, 0.5f));
	}
	if(ceiling) {
		// construct ceiling
		boxes.push_back(box3(vec3(-0.5f*w - W, h, -0.5f*d - W), vec3(0.5f*w + W, h + W, 0.5f*d + W)));
		box_colors.push_back(rgb(0.5f, 0.5f, 0.8f));
	}
}

/// construct boxes for environment
void vr_cobotics::construct_environment(float s, float ew, float ed, float w, float d, float h) {
	std::default_random_engine generator;
	std::uniform_real_distribution<float> distribution(0, 1);
	unsigned n = unsigned(ew / s);
	unsigned m = unsigned(ed / s);
	float ox = 0.5f*float(n)*s;
	float oz = 0.5f*float(m)*s;
	for(unsigned i = 0; i < n; ++i) {
		float x = i * s - ox;
		for(unsigned j = 0; j < m; ++j) {
			float z = j * s - oz;
			if(fabsf(x) < 0.5f*w && fabsf(x + s) < 0.5f*w && fabsf(z) < 0.5f*d && fabsf(z + s) < 0.5f*d)
				continue;
			float h = 0.2f*(std::max(abs(x) - 0.5f*w, 0.0f) + std::max(abs(z) - 0.5f*d, 0.0f))*distribution(generator) + 0.1f;
			boxes.push_back(box3(vec3(x, 0.0f, z), vec3(x + s, h, z + s)));
			rgb color = cgv::media::color<float, cgv::media::HLS>(distribution(generator), 0.1f*distribution(generator) + 0.15f, 0.3f);
			box_colors.push_back(color);
			/*box_colors.push_back(
				rgb(0.3f*distribution(generator) + 0.3f,
					0.3f*distribution(generator) + 0.2f,
					0.2f*distribution(generator) + 0.1f));*/
		}
	}
}

/// construct boxes that can be moved around
void vr_cobotics::construct_movable_boxes(float tw, float td, float th, float tW, size_t nr) {
	/*
	vec3 extent(0.75f, 0.5f, 0.05f);
	movable_boxes.push_back(box3(-0.5f * extent, 0.5f * extent));
	movable_box_colors.push_back(rgb(0, 0, 0));
	movable_box_translations.push_back(vec3(0, 1.2f, 0));
	movable_box_rotations.push_back(quat(1, 0, 0, 0));
	*/
	std::default_random_engine generator;
	std::uniform_real_distribution<float> distribution(0, 1);
	std::uniform_real_distribution<float> signed_distribution(-1, 1);
	for(size_t i = 0; i < nr; ++i) {
		float x = distribution(generator);
		float y = distribution(generator);
		vec3 extent(distribution(generator), distribution(generator), distribution(generator));
		extent += 0.01f;
		extent *= std::min(tw, td)*0.1f;

		vec3 center(-0.5f*tw + x * tw, th + tW, -0.5f*td + y * td);
		movable_boxes.push_back(box3(-0.5f*extent, 0.5f*extent));
		movable_box_colors.push_back(rgb(distribution(generator), distribution(generator), distribution(generator)));
		movable_box_translations.push_back(center);
		quat rot(signed_distribution(generator), signed_distribution(generator), signed_distribution(generator), signed_distribution(generator));
		rot.normalize();
		movable_box_rotations.push_back(rot);
	}
}
/// construct a trash bin 
void vr_cobotics::construct_trash_bin(float cw, float cd, float ch, float cH, float x, float y, float z)
{
	rgb trash_bin_clr(0.8f, 0.7f, 0.7f);
	boxes.push_back(box3(vec3(-0.5f * cw + x, 0 + y, -0.5f * cd + z), vec3(0.5f * cw + x, ch + y, 0.5f * cd + z)));
	box_colors.push_back(trash_bin_clr);

	boxes.push_back(box3(vec3(-0.5f * cw + x, 0 + y, -0.5f * cd + z), vec3(-0.5f * cw - ch + x, cH + y, 0.5f * cd + z)));
	box_colors.push_back(trash_bin_clr);

	boxes.push_back(box3(vec3(0.5f * cw - ch + x, 0 + y, -0.5f * cd + z), vec3(0.5f * cw + x, cH + y, 0.5f * cd + z)));
	box_colors.push_back(trash_bin_clr);

	boxes.push_back(box3(vec3(-0.5f * cw + x, 0 + y, -0.5f * cd + z), vec3(0.5f * cw + x, cH + y, -0.5f * cd + ch + z)));
	box_colors.push_back(trash_bin_clr);

	boxes.push_back(box3(vec3(-0.5f * cw - ch + x, 0 + y, 0.5f * cd + ch + z), vec3(0.5f * cw + x, cH + y, 0.5f * cd + z)));
	box_colors.push_back(trash_bin_clr);
}
/// construct a scene with a table
void vr_cobotics::build_scene(float w, float d, float h, float W, float tw, float td, float th, float tW, float cw, float cd, float ch, float cH, float x, float y, float z)
{
	construct_room(w, d, h, W, false, false);
	construct_table(tw, td, th, tW);
	construct_environment(0.3f, 3 * w, 3 * d, w, d, h);
	//construct_environment(0.4f, 0.5f, 1u, w, d, h);
	construct_movable_boxes(tw, td, th, tW, 20);
	if (is_trashbin)
		construct_trash_bin(cw, cd, ch, cH, x, y ,z);
}

vr_cobotics::vr_cobotics() 
{
	frame_split = 0;
	extent_texcrd = vec2(0.5f, 0.5f);
	center_left  = vec2(0.5f,0.25f);
	center_right = vec2(0.5f,0.25f);
	seethrough_gamma = 0.33f;
	frame_width = frame_height = 0;
	background_distance = 2;
	background_extent = 2;
	undistorted = true;
	shared_texture = true;
	max_rectangle = false;
	nr_cameras = 0;
	camera_tex_id = -1;
	camera_aspect = 1;
	use_matrix = true;
	show_seethrough = false;
	is_trashbin = false;
	set_name("vr_cobotics");
	build_scene(5, 7, 3, 0.2f, 1.6f, 0.8f, 0.7f, 0.03f, 0.2f, 0.2f, 0.01f, 0.25f, 1.5f, 0.0f, 0.0f);
	vr_view_ptr = 0;
	ray_length = 2;
	last_kit_handle = 0;
	connect(cgv::gui::ref_vr_server().on_device_change, this, &vr_cobotics::on_device_change);

	mesh_scale = 0.0005f;
	mesh_location = dvec3(0, 0.85f, 0);
	mesh_orientation = dquat(1, 0, 0, 0);

	srs.radius = 0.005f;

	vr_events_stream = nullptr;
	box_trajectory_stream = nullptr;
	controller_trajectory_stream = nullptr;
	box_select_mode = false;
	trash_bin_select_mode = false;
	box_edit_mode = false;
	for (auto& a : grab_number) a = 0;
	new_box = box3(vec3(-0.05f), vec3(0.05f));
	new_box_color = rgb(88.f / 255.f, 24.f / 255.f, 69.f / 255.f);
	edit_box_selected_axis = 0;
	edit_box_step = 0.025f;
	edit_box_max_size = 0.2f;
	new_box_distance = 0.35f;

	label_outofdate = true;
	label_text = "Info Board";
	label_font_idx = 0;
	label_upright = true;
	label_face_type = cgv::media::font::FFA_BOLD;
	label_resolution = 256;
	label_size = 20.0f;
	label_color = rgb(1, 1, 1);

	cgv::media::font::enumerate_font_names(font_names);
	font_enum_decl = "enums='";
	for (unsigned i = 0; i < font_names.size(); ++i) {
		if (i>0)
			font_enum_decl += ";";
		std::string fn(font_names[i]);
		if (cgv::utils::to_lower(fn) == "calibri") {
			label_font_face = cgv::media::font::find_font(fn)->get_font_face(label_face_type);
			label_font_idx = i;
		}
		font_enum_decl += std::string(fn);
	}
	font_enum_decl += "'";
	state[0] = state[1] = state[2] = state[3] = IS_NONE;
}
	
void vr_cobotics::stream_help(std::ostream& os) {
	os << "vr_cobotics: no shortcuts defined" << std::endl;
}
	
void vr_cobotics::on_set(void* member_ptr)
{
	if (member_ptr == &label_face_type || member_ptr == &label_font_idx) {
		label_font_face = cgv::media::font::find_font(font_names[label_font_idx])->get_font_face(label_face_type);
		label_outofdate = true;
	}
	if ((member_ptr >= &label_color && member_ptr < &label_color + 1) ||
		member_ptr == &label_size || member_ptr == &label_text) {
		label_outofdate = true;
	}
	if (member_ptr == &log_vr_events && vr_events_stream) {
		if (!log_vr_events) { //close file
			vr_events_stream->close();
			box_trajectory_stream->close();
			controller_trajectory_stream->close();
			vr_events_stream = nullptr;
			box_trajectory_stream = nullptr;
			controller_trajectory_stream = nullptr;
			vr_events_record_path = "";
		}
		else { //start timer
			vrr_t_start = std::chrono::steady_clock::now();
		}
	}
	/*if (member_ptr == &listen_address) {

	}*/
	update_member(member_ptr);
	post_redraw();
}
	
bool vr_cobotics::handle(cgv::gui::event& e)
{
	// check if vr event flag is not set and don't process events in this case
	if ((e.get_flags() & cgv::gui::EF_VR) == 0)
		return false;

	// record controller events
	if (log_vr_events && vr_events_stream) {
		auto now = std::chrono::steady_clock::now();
		double t = time_stamp(vrr_t_start);
		//double t = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(now - vrr_t_start).count()) / 1e6;
		*vr_events_stream << t << " \""; //timestamp
		e.stream_out(*vr_events_stream); //raw data
		*vr_events_stream << "\"\n";
	}

	// check event id
	switch (e.get_kind()) {
	case cgv::gui::EID_KEY:
	{
		cgv::gui::vr_key_event& vrke = static_cast<cgv::gui::vr_key_event&>(e);
		if (vrke.get_action() != cgv::gui::KA_RELEASE) {
			switch (vrke.get_key()) {
			case vr::VR_LEFT_BUTTON1:
				std::cout << "button 1 of left controller pressed" << std::endl;
			case vr::VR_RIGHT_BUTTON1:
				return true;
			case vr::VR_LEFT_BUTTON2:
				std::cout << "button 2 of left controller pressed" << std::endl;
			case vr::VR_RIGHT_BUTTON2:
				return true;
			case vr::VR_LEFT_BUTTON0:
				std::cout << "button 0 of left controller pressed" << std::endl;
			case vr::VR_RIGHT_BUTTON0:
			{
				//set grabed box as template for new boxes
				int ci = vrke.get_controller_index();
				if (box_edit_mode && state[ci] == IS_OVER || state[ci] == IS_GRAB) {
					// iterate intersection points of current controller
					for (size_t i = 0; i < intersection_points.size(); ++i) {
						if (intersection_controller_indices[i] != ci)
							continue;
						// extract box index
						int bi = intersection_box_indices[i];
						int axis = edit_box_selected_axis;
						new_box = movable_boxes[bi];
						new_box_color = movable_box_colors[bi];
						label_outofdate = true;
						break; //get only the first box
					}
					post_redraw();
				}
				return true;
			}
			case vr::VR_LEFT_STICK_UP:
				std::cout << "touch pad of right controller pressed at up direction" << std::endl;
			case vr::VR_RIGHT_STICK_UP:
			{
				return true;
			}
			case vr::VR_LEFT_STICK_DOWN:
				std::cout << "touch pad of left controller pressed at down direction" << std::endl;
			case vr::VR_RIGHT_STICK_DOWN:
			{
				return true;
			}
			case vr::VR_LEFT_STICK_LEFT:
				std::cout << "touch pad of left controller pressed at left direction" << std::endl;
			case vr::VR_RIGHT_STICK_LEFT:
				return true;
			case vr::VR_RIGHT_STICK_RIGHT:
			case vr::VR_LEFT_STICK_RIGHT:
				return true;
			}
		}
		else {
			switch (vrke.get_key()) {
			case vr::VR_LEFT_BUTTON0:
				std::cout << "button 0 of left controller released" << std::endl;
			case vr::VR_RIGHT_BUTTON0:

				return true;
			}
		}
		break;
	}
	case cgv::gui::EID_THROTTLE:
	{
		cgv::gui::vr_throttle_event& vrte = static_cast<cgv::gui::vr_throttle_event&>(e);
		std::cout << "throttle " << vrte.get_throttle_index() << " of controller " << vrte.get_controller_index()
			<< " adjusted from " << vrte.get_last_value() << " to " << vrte.get_value() << std::endl;
		if (vrte.get_value() == 1.0f) {
			if (box_edit_mode) {
				//create a new box
				int ci = vrte.get_controller_index();
				if (state[ci] == IS_NONE || state[ci] == IS_OVER) {
					const vec3 up = vec3(0.f, 1.f, 0.f);
					vec3 origin, direction;
					vrte.get_state().controller[ci].put_ray(&origin(0), &direction(0));
					movable_box_colors.emplace_back(new_box_color);
					movable_box_rotations.emplace_back(quat(cross(up, direction), acos(dot(up, direction))));
					movable_box_translations.emplace_back(origin + direction * new_box_distance);
					movable_boxes.emplace_back(new_box);
					post_redraw();
				}
				else if (state[ci] == IS_GRAB) {
					edit_box_selected_axis = (edit_box_selected_axis + 1) % 3;
				}
				return true;
			}
			if (box_select_mode) {
				int ci = vrte.get_controller_index();
				if (state[ci] == IS_NONE || state[ci] == IS_OVER) {
					vec3 origin, direction;
					vrte.get_state().controller[ci].put_ray(&origin(0), &direction(0));
					compute_intersections(origin, direction, ci, ci == 0 ? rgb(1, 0, 0) : rgb(0, 0, 1));
					unsigned bi = 0;
					if (intersection_points.size() != 0) {
						for (size_t i = 0; i < intersection_points.size(); ++i) {
							if (intersection_controller_indices[i] != ci)
								continue;
							// extract box index
							bi = intersection_box_indices[i];
						}
					}
					send_selection(bi);
					std::cout << "send a message!" << std::endl;
				}
			}
			if (trash_bin_select_mode) {
				int ci = vrte.get_controller_index();
				if (state[ci] == IS_NONE || state[ci] == IS_OVER) {
					vec3 origin, direction;
					vrte.get_state().controller[ci].put_ray(&origin(0), &direction(0));
					compute_intersections(origin, direction, ci, ci == 0 ? rgb(1, 0, 0) : rgb(0, 0, 1));
					unsigned bi = 0;
					if (intersection_points.size() != 0) {
						for (size_t i = 0; i < intersection_points.size(); ++i) {
							if (intersection_controller_indices[i] != ci)
								continue;
							// extract box index
							bi = intersection_box_indices[i];
						}
					}
					send_selection(bi);
				}
				std::cout << "send a message about trash bin!" << std::endl;
			}
		}
		return true;
	}
	case cgv::gui::EID_STICK:
	{
		cgv::gui::vr_stick_event& vrse = static_cast<cgv::gui::vr_stick_event&>(e);
		switch (vrse.get_action()) {
		case cgv::gui::SA_TOUCH:
			if (state[vrse.get_controller_index()] == IS_OVER)
				state[vrse.get_controller_index()] = IS_GRAB;
			break;
		case cgv::gui::SA_RELEASE:
			if (state[vrse.get_controller_index()] == IS_GRAB) {
				state[vrse.get_controller_index()] = IS_OVER;
				++grab_number[vrse.get_controller_index()];
			}
			break;
		case cgv::gui::SA_PRESS:
		{
			int ci = vrse.get_controller_index();
			if (box_edit_mode && vrse.get_x() < -0.5f && state[ci] == IS_GRAB) {
				//change box extents
				change_box_extents(AXIS_X, ci);
				post_redraw();
				return true;
			}
			if (box_edit_mode && vrse.get_y() < -0.5f && state[ci] == IS_GRAB) {
				//change box extents
				change_box_extents(AXIS_Y, ci);
				post_redraw();
				return true;
			}
			if (box_edit_mode && vrse.get_x() > 0.5f && state[ci] == IS_GRAB) {
				//change box extents
				change_box_extents(AXIS_Z, ci);
				post_redraw();
				return true;
			}
			if (box_edit_mode && vrse.get_y() > 0.5f && state[ci] == IS_GRAB) {
				//TODO delete box
				
				for (size_t i = 0; i < intersection_points.size(); ++i) {
					if (intersection_controller_indices[i] != ci)
						continue;
					// extract box index
					delete_box(intersection_box_indices[i]);
					break; //delete only the first box
				}
				// clear intersections of current controller 
				size_t i = clear_intersections(ci);
				// compute intersections
				vec3 origin, direction;
				vrse.get_state().controller[ci].put_ray(&origin(0), &direction(0));
				compute_intersections(origin, direction, ci, ci == 0 ? rgb(1, 0, 0) : rgb(0, 0, 1));
				label_outofdate = true;


				// update state based on whether we have found at least 
				// one intersection with controller ray
				if (intersection_points.size() == i)
					state[ci] = IS_NONE;
				else
					if (state[ci] == IS_NONE)
						state[ci] = IS_OVER;

				post_redraw();
				return true;
			}
		}
		case cgv::gui::SA_UNPRESS:
			std::cout << "stick " << vrse.get_stick_index()
				<< " of controller " << vrse.get_controller_index()
				<< " " << cgv::gui::get_stick_action_string(vrse.get_action())
				<< " at " << vrse.get_x() << ", " << vrse.get_y() << std::endl;
			return true;
		case cgv::gui::SA_MOVE:
		case cgv::gui::SA_DRAG:
			return true;
			std::cout << "stick " << vrse.get_stick_index()
				<< " of controller " << vrse.get_controller_index()
				<< " " << cgv::gui::get_stick_action_string(vrse.get_action())
				<< " from " << vrse.get_last_x() << ", " << vrse.get_last_y()
				<< " to " << vrse.get_x() << ", " << vrse.get_y() << std::endl;
			return true;
		}
		return true;
	}
	case cgv::gui::EID_POSE:
		cgv::gui::vr_pose_event& vrpe = static_cast<cgv::gui::vr_pose_event&>(e);
		// check for controller pose events
		int ci = vrpe.get_trackable_index();
		if (ci != -1) {
			if (state[ci] == IS_GRAB) {
				// in grab mode apply relative transformation to grabbed boxes

				// get previous and current controller position
				vec3 last_pos = vrpe.get_last_position();
				vec3 pos = vrpe.get_position();
				// get rotation from previous to current orientation
				// this is the current orientation matrix times the
				// inverse (or transpose) of last orientation matrix:
				// vrpe.get_orientation()*transpose(vrpe.get_last_orientation())
				mat3 rotation = vrpe.get_rotation_matrix();
				// iterate intersection points of current controller
				for (size_t i = 0; i < intersection_points.size(); ++i) {
					if (intersection_controller_indices[i] != ci)
						continue;
					// extract box index
					unsigned bi = intersection_box_indices[i];
					// update translation with position change and rotation
					movable_box_translations[bi] = 
						rotation * (movable_box_translations[bi] - last_pos) + pos;
					// update orientation with rotation, note that quaternions
					// need to be multiplied in oposite order. In case of matrices
					// one would write box_orientation_matrix *= rotation
					movable_box_rotations[bi] = quat(rotation) * movable_box_rotations[bi];
					// update intersection points
					intersection_points[i] = rotation * (intersection_points[i] - last_pos) + pos;

					if (log_vr_events && box_trajectory_stream) {
						//<box-index> <controller-id> <grab-id> <box-translation> <box-rotation>
						*box_trajectory_stream
							<< ci << " "
							<< grab_number[ci] << " "
							<< time_stamp(vrr_t_start) << " "
							<< bi << " "
							<< movable_box_translations[bi] << " "
							<< movable_box_rotations[bi] << '\n';
					}
				}
			}
			else {// not grab
				size_t i = clear_intersections(ci);

				// compute intersections
				vec3 origin, direction;
				vrpe.get_state().controller[ci].put_ray(&origin(0), &direction(0));
				compute_intersections(origin, direction, ci, ci == 0 ? rgb(1, 0, 0) : rgb(0, 0, 1));
				label_outofdate = true;


				// update state based on whether we have found at least 
				// one intersection with controller ray
				if (intersection_points.size() == i)
					state[ci] = IS_NONE;
				else
					if (state[ci] == IS_NONE)
						state[ci] = IS_OVER;
			}
			//log controller trajectory
			if (log_vr_events && controller_trajectory_stream) {
				mat3 rotation = vrpe.get_rotation_matrix();
				*controller_trajectory_stream
					<< ci << " "
					<< ((state[ci] == IS_GRAB) ? grab_number[ci] : -1) << " "
					<< time_stamp(vrr_t_start) << " "
					<< vrpe.get_position() << " "
					<< quat(vrpe.get_orientation()) << '\n';
			}
			post_redraw();
		}
		return true;
	}
	return false;
}

bool vr_cobotics::init(cgv::render::context& ctx)
{
	if (!cgv::utils::has_option("NO_OPENVR"))
		ctx.set_gamma(1.0f);

	if (!seethrough.build_program(ctx, "seethrough.glpr"))
		cgv::gui::message("could not build seethrough program");
	
	cgv::media::mesh::simple_mesh<> M;
#ifdef _DEBUG
	if (M.read("D:/data/surface/meshes/obj/Max-Planck_lowres.obj")) {
#else
	if (M.read("D:/data/surface/meshes/obj/Max-Planck_highres.obj")) {
#endif
		MI.construct(ctx, M);
		MI.bind(ctx, ctx.ref_surface_shader_program(true), true);
	}

	cgv::gui::connect_vr_server(true);

	auto view_ptr = find_view_as_node();
	if (view_ptr) {
		view_ptr->set_eye_keep_view_angle(dvec3(0, 4, -4));
		// if the view points to a vr_view_interactor
		vr_view_ptr = dynamic_cast<vr_view_interactor*>(view_ptr);
		if (vr_view_ptr) {
			// configure vr event processing
			vr_view_ptr->set_event_type_flags(
				cgv::gui::VREventTypeFlags(
					cgv::gui::VRE_KEY +
					cgv::gui::VRE_THROTTLE +
					cgv::gui::VRE_STICK +
					cgv::gui::VRE_STICK_KEY +
					cgv::gui::VRE_POSE
				));
			vr_view_ptr->enable_vr_event_debugging(false);
			// configure vr rendering
			vr_view_ptr->draw_action_zone(false);
			vr_view_ptr->draw_vr_kits(true);
			vr_view_ptr->enable_blit_vr_views(true);
			vr_view_ptr->set_blit_vr_view_width(200);
		}
	}

	cgv::render::ref_box_renderer(ctx, 1);
	cgv::render::ref_sphere_renderer(ctx, 1);
	cgv::render::ref_rounded_cone_renderer(ctx, 1);
	return true;
}

void vr_cobotics::clear(cgv::render::context& ctx)
{
	cgv::render::ref_box_renderer(ctx, -1);
	cgv::render::ref_sphere_renderer(ctx, -1);
	cgv::render::ref_rounded_cone_renderer(ctx, -1);
}

void vr_cobotics::init_frame(cgv::render::context& ctx)
{
	if (label_fbo.get_width() != label_resolution) {
		label_tex.destruct(ctx);
		label_fbo.destruct(ctx);
	}
	if (!label_fbo.is_created()) {
		label_tex.create(ctx, cgv::render::TT_2D, label_resolution, label_resolution);
		label_fbo.create(ctx, label_resolution, label_resolution);
		label_tex.set_min_filter(cgv::render::TF_LINEAR_MIPMAP_LINEAR);
		label_tex.set_mag_filter(cgv::render::TF_LINEAR);
		label_fbo.attach(ctx, label_tex);
		label_outofdate = true;
	}
	if (label_outofdate && label_fbo.is_complete(ctx)) {
		glPushAttrib(GL_COLOR_BUFFER_BIT);
		label_fbo.enable(ctx);
		label_fbo.push_viewport(ctx);
		ctx.push_pixel_coords();
			glClearColor(0.5f,0.5f,0.5f,1.0f);
			glClear(GL_COLOR_BUFFER_BIT);

			glColor4f(label_color[0], label_color[1], label_color[2], 1);
			ctx.set_cursor(20, (int)ceil(label_size) + 20);
			ctx.enable_font_face(label_font_face, label_size);
			ctx.output_stream() << label_text << "\n";
			ctx.output_stream().flush(); // make sure to flush the stream before change of font size or font face

			ctx.enable_font_face(label_font_face, 0.7f*label_size);
			
			if (box_edit_mode) {
				const char axis[] = "XYZ";
				ctx.output_stream() << "new box[Trigger] \nextent=" << new_box.get_extent() << "\ncolor=" << new_box_color << '\n';
				ctx.output_stream() << "set as box template[Grip Button]\n";
				ctx.output_stream() << "delete box[trackpad up]\n";
				ctx.output_stream() << "change box extents[trackpad left,right,down]\n";
			}
			
			for (size_t i = 0; i < intersection_points.size(); ++i) {
				ctx.output_stream()
					<< "box " << intersection_box_indices[i]
					<< " at (" << intersection_points[i]
					<< ") with controller " << intersection_controller_indices[i] << "\n";
			}
			ctx.output_stream().flush();

		ctx.pop_pixel_coords();
		label_fbo.pop_viewport(ctx);
		label_fbo.disable(ctx);
		glPopAttrib();
		label_outofdate = false;

		label_tex.generate_mipmaps(ctx);
	}
	if (vr_view_ptr && vr_view_ptr->get_rendered_vr_kit() != 0 && vr_view_ptr->get_rendered_eye() == 0 && vr_view_ptr->get_rendered_vr_kit() == vr_view_ptr->get_current_vr_kit()) {
		vr::vr_kit* kit_ptr = vr_view_ptr->get_current_vr_kit();
		if (kit_ptr) {
			vr::vr_camera* camera_ptr = kit_ptr->get_camera();
			if (camera_ptr && camera_ptr->get_state() == vr::CS_STARTED) {
				uint32_t width = frame_width, height = frame_height, split = frame_split;
				if (shared_texture) {
					box2 tex_range;
					if (camera_ptr->get_gl_texture_id(camera_tex_id, width, height, undistorted, &tex_range.ref_min_pnt()(0))) {
						camera_aspect = (float)width / height;
						split = camera_ptr->get_frame_split();
						switch (split) {
						case vr::CFS_VERTICAL:
							camera_aspect *= 2;
							break;
						case vr::CFS_HORIZONTAL:
							camera_aspect *= 0.5f;
							break;
						}
					}
					else
						camera_tex_id = -1;
				}
				else {
					std::vector<uint8_t> frame_data;
					if (camera_ptr->get_frame(frame_data, width, height, undistorted, max_rectangle)) {
						camera_aspect = (float)width / height;
						split = camera_ptr->get_frame_split();
						switch (split) {
						case vr::CFS_VERTICAL:
							camera_aspect *= 2;
							break;
						case vr::CFS_HORIZONTAL:
							camera_aspect *= 0.5f;
							break;
						}
						cgv::data::data_format df(width, height, cgv::type::info::TI_UINT8, cgv::data::CF_RGBA);
						cgv::data::data_view dv(&df, frame_data.data());
						if (camera_tex.is_created()) {
							if (camera_tex.get_width() != width || camera_tex.get_height() != height)
								camera_tex.destruct(ctx);
							else
								camera_tex.replace(ctx, 0, 0, dv);
						}
						if (!camera_tex.is_created())
							camera_tex.create(ctx, dv);
					}
					else if (camera_ptr->has_error())
						cgv::gui::message(camera_ptr->get_last_error());
				}
				if (frame_width != width || frame_height != height) {
					frame_width = width;
					frame_height = height;

					center_left(0) = camera_centers[2](0) / frame_width;
					center_left(1) = camera_centers[2](1) / frame_height;
					center_right(0) = camera_centers[3](0) / frame_width;
					center_right(1) = camera_centers[3](1) / frame_height;

					update_member(&frame_width);
					update_member(&frame_height);
					update_member(&center_left(0));
					update_member(&center_left(1));
					update_member(&center_right(0));
					update_member(&center_right(1));
				}
				if (split != frame_split) {
					frame_split = split;
					update_member(&frame_split);
				}
			}
		}
	}
}

void vr_cobotics::draw(cgv::render::context& ctx)
{
	if (MI.is_constructed()) {
		dmat4 R;
		mesh_orientation.put_homogeneous_matrix(R);
		ctx.push_modelview_matrix();
		ctx.mul_modelview_matrix(
			cgv::math::translate4<double>(mesh_location)*
			cgv::math::scale4<double>(mesh_scale, mesh_scale, mesh_scale) *
			R);
		MI.draw_all(ctx);
		ctx.pop_modelview_matrix();
	}
	if (vr_view_ptr) {
		if ((!shared_texture && camera_tex.is_created()) || (shared_texture && camera_tex_id != -1)) {
			if (vr_view_ptr->get_rendered_vr_kit() != 0 && vr_view_ptr->get_rendered_vr_kit() == vr_view_ptr->get_current_vr_kit()) {
				int eye = vr_view_ptr->get_rendered_eye();

				// compute billboard
				dvec3 vd = vr_view_ptr->get_view_dir_of_kit();
				dvec3 y = vr_view_ptr->get_view_up_dir_of_kit();
				dvec3 x = normalize(cross(vd, y));
				y = normalize(cross(x, vd));
				x *= camera_aspect * background_extent * background_distance;
				y *= background_extent * background_distance;
				vd *= background_distance;
				dvec3 eye_pos = vr_view_ptr->get_eye_of_kit(eye);
				std::vector<vec3> P;
				std::vector<vec2> T;
				P.push_back(eye_pos + vd - x - y);
				P.push_back(eye_pos + vd + x - y);
				P.push_back(eye_pos + vd - x + y);
				P.push_back(eye_pos + vd + x + y);
				double v_offset = 0.5 * (1 - eye);
				T.push_back(dvec2(0.0, 0.5 + v_offset));
				T.push_back(dvec2(1.0, 0.5 + v_offset));
				T.push_back(dvec2(0.0, v_offset));
				T.push_back(dvec2(1.0, v_offset));

				cgv::render::shader_program& prog = seethrough;
				cgv::render::attribute_array_binding::set_global_attribute_array(ctx, prog.get_position_index(), P);
				cgv::render::attribute_array_binding::set_global_attribute_array(ctx, prog.get_texcoord_index(), T);
				cgv::render::attribute_array_binding::enable_global_array(ctx, prog.get_position_index());
				cgv::render::attribute_array_binding::enable_global_array(ctx, prog.get_texcoord_index());

				GLint active_texture, texture_binding;
				if (shared_texture) {
					glGetIntegerv(GL_ACTIVE_TEXTURE, &active_texture);
					glGetIntegerv(GL_TEXTURE_BINDING_2D, &texture_binding);
					glActiveTexture(GL_TEXTURE0);
					glBindTexture(GL_TEXTURE_2D, camera_tex_id);
				}
				else
					camera_tex.enable(ctx, 0);
				prog.set_uniform(ctx, "texture", 0);
				prog.set_uniform(ctx, "seethrough_gamma", seethrough_gamma);
				prog.set_uniform(ctx, "use_matrix", use_matrix);

				// use of convenience function
				vr::configure_seethrough_shader_program(ctx, prog, frame_width, frame_height,
					vr_view_ptr->get_current_vr_kit(), *vr_view_ptr->get_current_vr_state(),
					0.01f, 2 * background_distance, eye, undistorted);

				/* equivalent detailed code relies on more knowledge on program parameters
				mat4 TM = vr::get_texture_transform(vr_view_ptr->get_current_vr_kit(), *vr_view_ptr->get_current_vr_state(), 0.01f, 2 * background_distance, eye, undistorted);
				prog.set_uniform(ctx, "texture_matrix", TM);

				prog.set_uniform(ctx, "extent_texcrd", extent_texcrd);
				prog.set_uniform(ctx, "frame_split", frame_split);
				prog.set_uniform(ctx, "center_left", center_left);
				prog.set_uniform(ctx, "center_right", center_right);
				prog.set_uniform(ctx, "eye", eye);
				*/
				prog.enable(ctx);
				ctx.set_color(rgba(1, 1, 1, 1));

				glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);


				prog.disable(ctx);

				if (shared_texture) {
					glActiveTexture(active_texture);
					glBindTexture(GL_TEXTURE_2D, texture_binding);
				}
				else
					camera_tex.disable(ctx);

				cgv::render::attribute_array_binding::disable_global_array(ctx, prog.get_position_index());
				cgv::render::attribute_array_binding::disable_global_array(ctx, prog.get_texcoord_index());
			}
		}
		if (vr_view_ptr) {
			std::vector<vec3> P;
			std::vector<float> R;
			std::vector<rgb> C;
			const vr::vr_kit_state* state_ptr = vr_view_ptr->get_current_vr_state();
			if (state_ptr) {
				for (int ci = 0; ci < 4; ++ci) if (state_ptr->controller[ci].status == vr::VRS_TRACKED) {
					vec3 ray_origin, ray_direction;
					state_ptr->controller[ci].put_ray(&ray_origin(0), &ray_direction(0));
					P.push_back(ray_origin);
					R.push_back(0.002f);
					P.push_back(ray_origin + ray_length * ray_direction);
					R.push_back(0.003f);
					rgb c(float(1 - ci), 0.5f * (int)state[ci], float(ci));
					C.push_back(c);
					C.push_back(c);
				}
			}
			if (P.size() > 0) {
				auto& cr = cgv::render::ref_rounded_cone_renderer(ctx);
				cr.set_render_style(cone_style);
				//cr.set_eye_position(vr_view_ptr->get_eye_of_kit());
				cr.set_position_array(ctx, P);
				cr.set_color_array(ctx, C);
				cr.set_radius_array(ctx, R);
				if (!cr.render(ctx, 0, P.size())) {
					cgv::render::shader_program& prog = ctx.ref_default_shader_program();
					int pi = prog.get_position_index();
					int ci = prog.get_color_index();
					cgv::render::attribute_array_binding::set_global_attribute_array(ctx, pi, P);
					cgv::render::attribute_array_binding::enable_global_array(ctx, pi);
					cgv::render::attribute_array_binding::set_global_attribute_array(ctx, ci, C);
					cgv::render::attribute_array_binding::enable_global_array(ctx, ci);
					glLineWidth(3);
					prog.enable(ctx);
					glDrawArrays(GL_LINES, 0, (GLsizei)P.size());
					prog.disable(ctx);
					cgv::render::attribute_array_binding::disable_global_array(ctx, pi);
					cgv::render::attribute_array_binding::disable_global_array(ctx, ci);
					glLineWidth(1);
				}
			}
		}
	}
	cgv::render::box_renderer& renderer = cgv::render::ref_box_renderer(ctx);

	// draw wireframe boxes
	if (frame_boxes.size() > 0) { //pointing the renderer to uninitialized arrays causes assertions in debug mode
		renderer.set_render_style(wire_frame_style);
		renderer.set_box_array(ctx, frame_boxes);
		renderer.set_color_array(ctx, frame_box_colors);
		renderer.set_translation_array(ctx, frame_box_translations);
		renderer.set_rotation_array(ctx, frame_box_rotations);
		if (renderer.validate_and_enable(ctx)) {
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			renderer.draw(ctx, 0, frame_boxes.size());
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}
		renderer.disable(ctx);
	}

	// draw dynamic boxes 
	renderer.set_render_style(movable_style);
	renderer.set_box_array(ctx, movable_boxes);
	renderer.set_color_array(ctx, movable_box_colors);
	renderer.set_translation_array(ctx, movable_box_translations);
	renderer.set_rotation_array(ctx, movable_box_rotations);
	if (renderer.validate_and_enable(ctx)) {
		if (show_seethrough) {
			glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
			renderer.draw(ctx, 0, 3);
			glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
			renderer.draw(ctx, 3, movable_boxes.size() - 3);
		}
		else
			renderer.draw(ctx, 0, movable_boxes.size());
	}
	renderer.disable(ctx);

	// draw static boxes
	renderer.set_render_style(style);
	renderer.set_box_array(ctx, boxes);
	renderer.set_color_array(ctx, box_colors);
	renderer.render(ctx, 0, boxes.size());


	// draw intersection points
	if (!intersection_points.empty()) {
		auto& sr = cgv::render::ref_sphere_renderer(ctx);
		sr.set_position_array(ctx, intersection_points);
		sr.set_color_array(ctx, intersection_colors);
		sr.set_render_style(srs);
		sr.render(ctx, 0, intersection_points.size());
	}

	// draw label
	if (vr_view_ptr && label_tex.is_created()) {
		cgv::render::shader_program& prog = ctx.ref_default_shader_program(true);
		int pi = prog.get_position_index();
		int ti = prog.get_texcoord_index();
		vec3 p(0, 1.5f, 0);
		vec3 y = label_upright ? vec3(0, 1.0f, 0) : normalize(vr_view_ptr->get_view_up_dir_of_kit());
		vec3 x = normalize(cross(vec3(vr_view_ptr->get_view_dir_of_kit()), y));
		float w = 0.5f, h = 0.5f;
		std::vector<vec3> P;
		std::vector<vec2> T;
		P.push_back(p - 0.5f * w * x - 0.5f * h * y); T.push_back(vec2(0.0f, 0.0f));
		P.push_back(p + 0.5f * w * x - 0.5f * h * y); T.push_back(vec2(1.0f, 0.0f));
		P.push_back(p - 0.5f * w * x + 0.5f * h * y); T.push_back(vec2(0.0f, 1.0f));
		P.push_back(p + 0.5f * w * x + 0.5f * h * y); T.push_back(vec2(1.0f, 1.0f));
		cgv::render::attribute_array_binding::set_global_attribute_array(ctx, pi, P);
		cgv::render::attribute_array_binding::enable_global_array(ctx, pi);
		cgv::render::attribute_array_binding::set_global_attribute_array(ctx, ti, T);
		cgv::render::attribute_array_binding::enable_global_array(ctx, ti);
		prog.enable(ctx);
		label_tex.enable(ctx);
		ctx.set_color(rgb(1, 1, 1));
		glDrawArrays(GL_TRIANGLE_STRIP, 0, (GLsizei)P.size());
		label_tex.disable(ctx);
		prog.disable(ctx);
		cgv::render::attribute_array_binding::disable_global_array(ctx, pi);
		cgv::render::attribute_array_binding::disable_global_array(ctx, ti);
	}
}

void vr_cobotics::finish_draw(cgv::render::context& ctx)
{
	return;
	if ((!shared_texture && camera_tex.is_created()) || (shared_texture && camera_tex_id != -1)) {
		cgv::render::shader_program& prog = ctx.ref_default_shader_program(true);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		GLint active_texture, texture_binding;
		if (shared_texture) {
			glGetIntegerv(GL_ACTIVE_TEXTURE, &active_texture);
			glGetIntegerv(GL_TEXTURE_BINDING_2D, &texture_binding);
			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, camera_tex_id);
		}
		else
			camera_tex.enable(ctx, 0);

		prog.set_uniform(ctx, "texture", 0);
		ctx.push_modelview_matrix();
		ctx.mul_modelview_matrix(cgv::math::translate4<double>(0, 3, 0));
		prog.enable(ctx);
		ctx.set_color(rgba(1, 1, 1, 0.8f));
		ctx.tesselate_unit_square();
		prog.disable(ctx);
		if (shared_texture) {
			glActiveTexture(active_texture);
			glBindTexture(GL_TEXTURE_2D, texture_binding);
		}
		else
			camera_tex.disable(ctx);
		ctx.pop_modelview_matrix();
		glDisable(GL_BLEND);
	}
}

					//0.2f*distribution(generator)+0.1f));
void vr_cobotics::create_gui() {
	add_decorator("vr_cobotics", "heading", "level=2");
	add_member_control(this, "mesh_scale", mesh_scale, "value_slider", "min=0.1;max=10;log=true;ticks=true");
	add_gui("mesh_location", mesh_location, "vector", "options='min=-3;max=3;ticks=true");
	add_gui("mesh_orientation", static_cast<dvec4&>(mesh_orientation), "direction", "options='min=-1;max=1;ticks=true");
	add_member_control(this, "ray_length", ray_length, "value_slider", "min=0.1;max=10;log=true;ticks=true");
	add_member_control(this, "show_seethrough", show_seethrough, "check");
	if(last_kit_handle) {
		add_decorator("cameras", "heading", "level=3");
		add_view("nr", nr_cameras);
		if(nr_cameras > 0) {
			connect_copy(add_button("start")->click, cgv::signal::rebind(this, &vr_cobotics::start_camera));
			connect_copy(add_button("stop")->click, cgv::signal::rebind(this, &vr_cobotics::stop_camera));
			add_view("frame_width", frame_width, "", "w=20", "  ");
			add_view("height", frame_height, "", "w=20", "  ");
			add_view("split", frame_split, "", "w=50");
			add_member_control(this, "undistorted", undistorted, "check");
			add_member_control(this, "shared_texture", shared_texture, "check");
			add_member_control(this, "max_rectangle", max_rectangle, "check");
			add_member_control(this, "use_matrix", use_matrix, "check");
			add_member_control(this, "gamma", seethrough_gamma, "value_slider", "min=0.1;max=10;log=true;ticks=true");
			add_member_control(this, "extent_x", extent_texcrd[0], "value_slider", "min=0.2;max=2;ticks=true");
			add_member_control(this, "extent_y", extent_texcrd[1], "value_slider", "min=0.2;max=2;ticks=true");
			add_member_control(this, "center_left_x", center_left[0], "value_slider", "min=0.2;max=0.8;ticks=true");
			add_member_control(this, "center_left_y", center_left[1], "value_slider", "min=0.2;max=0.8;ticks=true");
			add_member_control(this, "center_right_x", center_right[0], "value_slider", "min=0.2;max=0.8;ticks=true");
			add_member_control(this, "center_right_y", center_right[1], "value_slider", "min=0.2;max=0.8;ticks=true");
			add_member_control(this, "background_distance", background_distance, "value_slider", "min=0.1;max=10;log=true;ticks=true");
			add_member_control(this, "background_extent", background_extent, "value_slider", "min=0.01;max=10;log=true;ticks=true");
		}
		vr::vr_kit* kit_ptr = vr::get_vr_kit(last_kit_handle);
		const std::vector<std::pair<int, int> >* t_and_s_ptr = 0;
		if(kit_ptr)
			t_and_s_ptr = &kit_ptr->get_controller_throttles_and_sticks(0);
		add_decorator("deadzone and precisions", "heading", "level=3");
		int ti = 0;
		int si = 0;
		for(unsigned i = 0; i < left_deadzone_and_precision.size(); ++i) {
			std::string prefix = std::string("unknown[") + cgv::utils::to_string(i) + "]";
			if(t_and_s_ptr) {
				if(t_and_s_ptr->at(i).second == -1)
					prefix = std::string("throttle[") + cgv::utils::to_string(ti++) + "]";
				else
					prefix = std::string("stick[") + cgv::utils::to_string(si++) + "]";
			}
			add_member_control(this, prefix + ".deadzone", left_deadzone_and_precision[i].first, "value_slider", "min=0;max=1;ticks=true;log=true");
			add_member_control(this, prefix + ".precision", left_deadzone_and_precision[i].second, "value_slider", "min=0;max=1;ticks=true;log=true");
		}
	}
	if (begin_tree_node("NNG events", 1.f)) {
		align("\a");
		connect_copy(add_button("receive")->click, rebind(this, &vr_cobotics::on_send_movable_boxes_id_cb));
		add_member_control(this, "select box", box_select_mode, "toggle");
		add_member_control(this, "select trash bin", trash_bin_select_mode, "toggle");
		align("\b");
		end_tree_node(style);
	}
	if (begin_tree_node("movable boxes", box_edit_mode)) {
		align("\a");
		connect_copy(add_button("clear boxes")->click, rebind(this, &vr_cobotics::clear_movable_boxes));
		connect_copy(add_button("save boxes")->click, rebind(this, &vr_cobotics::on_save_movable_boxes_cb));
		connect_copy(add_button("load boxes")->click, rebind(this, &vr_cobotics::on_load_movable_boxes_cb));
		connect_copy(add_button("load target")->click, rebind(this, &vr_cobotics::on_load_wireframe_boxes_cb));
		add_member_control(this, "allow editing/creating boxes", box_edit_mode, "toggle");
		add_member_control(this, "max. box size", edit_box_max_size, "value_slider", "min=0;max=1;ticks=true");
		add_member_control(this, "edit step ize", edit_box_step, "value_slider", "min=0;max=1;ticks=true");
		end_tree_node(box_edit_mode);
		align("\b");
		end_tree_node(style);
	}
	if (begin_tree_node("VR events", log_vr_events)) {
		align("\a");
		add_member_control(this, "log vr events", log_vr_events, "toggle");
		connect_copy(add_button("select protocol file")->click, rebind(this, &vr_cobotics::on_set_vr_event_streaming_file));
		add_view("protocol file", vr_events_record_path);
		end_tree_node(log_vr_events);
		align("\b");
		end_tree_node(style);
	}
	if (begin_tree_node("box style", style)) {
		align("\a");
		add_gui("box style", style);
		align("\b");
		end_tree_node(style);
	}
	if (begin_tree_node("cone style", cone_style)) {
		align("\a");
		add_gui("cone style", cone_style);
		align("\b");
		end_tree_node(cone_style);
	}
	if(begin_tree_node("movable box style", movable_style)) {
		align("\a");
		add_gui("movable box style", movable_style);
		align("\b");
		end_tree_node(movable_style);
	}
	if(begin_tree_node("intersections", srs)) {
		align("\a");
		add_gui("sphere style", srs);
		align("\b");
		end_tree_node(srs);
	}
	if(begin_tree_node("mesh", mesh_scale)) {
		align("\a");
		add_member_control(this, "scale", mesh_scale, "value_slider", "min=0.0001;step=0.0000001;max=100;log=true;ticks=true");
		add_gui("location", mesh_location, "", "main_label='';long_label=true;gui_type='value_slider';options='min=-2;max=2;step=0.001;ticks=true'");
		add_gui("orientation", static_cast<dvec4&>(mesh_orientation), "direction", "main_label='';long_label=true;gui_type='value_slider';options='min=-1;max=1;step=0.001;ticks=true'");
		align("\b");
		end_tree_node(mesh_scale);
	}

	if(begin_tree_node("label", label_size)) {
		align("\a");
		add_member_control(this, "text", label_text);
		add_member_control(this, "upright", label_upright);
		add_member_control(this, "font", (cgv::type::DummyEnum&)label_font_idx, "dropdown", font_enum_decl);
		add_member_control(this, "face", (cgv::type::DummyEnum&)label_face_type, "dropdown", "enums='regular,bold,italics,bold+italics'");
		add_member_control(this, "size", label_size, "value_slider", "min=8;max=64;ticks=true");
		add_member_control(this, "color", label_color);
		add_member_control(this, "resolution", (cgv::type::DummyEnum&)label_resolution, "dropdown", "enums='256=256,512=512,1024=1024,2048=2048'");
		align("\b");
		end_tree_node(label_size);
	}
}

bool vr_cobotics::self_reflect(cgv::reflect::reflection_handler & rh)
{
	return	rh.reflect_member("vr_events_record_path", vr_events_record_path) &&
		    rh.reflect_member("select a box", box_select_mode) &&
			rh.reflect_member("select a trash bin", trash_bin_select_mode) &&
			rh.reflect_member("box_edit_mode", box_edit_mode);
}

bool vr_cobotics::save_boxes(const std::string fn, const std::vector<box3>& boxes, const std::vector<rgb>& box_colors, const std::vector<vec3>& box_translations, const std::vector<quat>& box_rotations)
{
	std::stringstream data;


	if (boxes.size() != box_colors.size() || boxes.size() != box_translations.size() || boxes.size() != box_rotations.size()) {
		std::cerr << "vr_cobotics::save_boxes: passed vectors have different sizes!";
		return false;
	}

	for (size_t i = 0; i < movable_boxes.size(); ++i) {
		//format: BOX <box.min_p> <box.max_p> <trans> <rot> <col>
		const vec3& box_translation = box_translations[i];
		const quat& box_rotation = box_rotations[i];
		const rgb& box_color = box_colors[i];
		data << "BOX "
			<< boxes[i].get_min_pnt() << " "
			<< boxes[i].get_max_pnt() << " "
			<< box_translation << " "
			<< box_rotation << " "
			<< box_color << " "
			<< '\n';
	}
	std::string s = data.str();
	if (!cgv::utils::file::write(fn, s.data(), s.size())) {
		std::cerr << "vr_cobotics::save_boxes: failed writing data to file: " << fn;
	}
	return true;
}

bool vr_cobotics::load_boxes(const std::string fn, std::vector<box3>& boxes, std::vector<rgb>& box_colors, std::vector<vec3>& box_translations, std::vector<quat>& box_rotations)
{
	std::string data;
	if (!cgv::utils::file::read(fn, data)) {
		std::cerr << "vr_cobotics::load_boxes: failed reading data from file: " << fn << '\n';
		return false;
	}
	std::istringstream f(data);
	std::string line;

	while (!f.eof()) {
		std::getline(f, line); 	//read a line
		std::istringstream l(line);
		std::string sym;
		
		int limit=1;
		bool valid = true;
		if (!l.eof()) {
			getline(l, sym, ' '); //get the first symbol determing the type
			if (sym == "BOX") { //in case of a box
				vec3 minp,maxp,trans;
				quat rot;
				rgb col;
				l >> minp;
				l >> maxp;
				l >> trans;
				l >> rot;
				l >> col;
				
				boxes.emplace_back(minp, maxp);
				box_translations.emplace_back(trans);
				box_rotations.emplace_back(rot);
				box_colors.emplace_back(col);
			}
		}
	}
	return true;
}

void vr_cobotics::resize_box(int box_index, vec3 extends)
{
	vec3 center = (boxes[box_index].get_max_pnt() + boxes[box_index].get_min_pnt())*0.5f;
	vec3 minp = center - 0.5f*extends;
	vec3 maxp = center + 0.5f*extends;
	boxes[box_index] = box3(minp, maxp);
}

void vr_cobotics::on_save_movable_boxes_cb()
{
	const std::string file_ending = ".vrboxes";
	std::string fn = cgv::gui::file_save_dialog("base file name", "Box configurations(vrboxes):*.vrboxes");
	if (!hasEnding(fn, file_ending)) {
		fn.append(file_ending);
	}
	if (fn.empty())
		return;
	
	save_boxes(fn, movable_boxes, movable_box_colors, movable_box_translations, movable_box_rotations);
}

void vr_cobotics::on_load_movable_boxes_cb()
{
	std::string fn = cgv::gui::file_open_dialog("base file name", "Box configurations(vrboxes):*.vrboxes");
	if (!cgv::utils::file::exists(fn)) {
		std::cerr << "vr_cobotics::on_load_movable_boxes_cb: file does not exist!\n";
		return;
	}
	clear_movable_boxes();
	if (!load_boxes(fn, movable_boxes, movable_box_colors, movable_box_translations, movable_box_rotations)) {
		std::cerr << "vr_cobotics::on_load_movable_boxes_cb: failed to parse file!\n";
		clear_movable_boxes(); //delete all boxes after a failure to reach a valid logical state
	}
}

void vr_cobotics::on_load_wireframe_boxes_cb()
{
	std::string fn = cgv::gui::file_open_dialog("base file name", "Box configurations(txt):*.txt");
	if (!cgv::utils::file::exists(fn)) {
		std::cerr << "vr_cobotics::on_load_movable_boxes_cb: file does not exist!\n";
		return;
	}
	clear_frame_boxes();
	if (!load_boxes(fn, frame_boxes, frame_box_colors, frame_box_translations, frame_box_rotations)) {
		std::cerr << "vr_cobotics::on_load_wireframe_boxes_cb: failed to parse file!\n";
		clear_frame_boxes(); //delete all boxes after a failure to reach a valid logical state
	}
}

void vr_cobotics::on_set_vr_event_streaming_file()
{
	std::string fn = cgv::gui::file_save_dialog("base file name", "File Prefix"); //VR-Conntroller Record
	if (fn.empty())
		return;
	vr_events_record_path = fn;

	vr_events_stream = std::make_shared<std::ofstream>(fn+".vrcr"); //VR Conntroller Record(.vrcr) :*.vrcr
	box_trajectory_stream = std::make_shared<std::ofstream>(fn + ".btrj"); //Block trajectory : *.btrj
	controller_trajectory_stream = std::make_shared<std::ofstream>(fn + ".ctrj"); //controller trajectory : *.ctrj
	if (!vr_events_stream->good()) {
		std::cerr << "vr_cobotics::on_set_vr_event_streaming_file: can't write file!\n";
		vr_events_stream = nullptr;
	}
	post_recreate_gui();
}

void vr_cobotics::on_send_movable_boxes_id_cb()
{
	sender();
}

void vr_cobotics::clear_movable_boxes()
{
	movable_boxes.clear();
	movable_box_translations.clear();
	movable_box_rotations.clear();
	movable_box_colors.clear();
}

void vr_cobotics::clear_frame_boxes()
{
	frame_boxes.clear();
	frame_box_translations.clear();
	frame_box_rotations.clear();
	frame_box_colors.clear();
}

Scene vr_cobotics::buildDummyScene()
{
	Scene scene;

	std::vector<std::string> names{ "object0", "object1", "object2", "object3", "bin" };

	std::vector<std::vector<float>> values{
		//   0-2: min_pos        3-5 max_pos      6-8 trans_XZY    9-12 quat.             13-15 rgb
		{-0.05, -0.05, -0.05, 0.05, 0.05, 0.05, 0.5,  0.75, -0.1, 1.0, 0, 0, 0, 0.221034, 0.913376, 0.968868},
		{-0.05, -0.05, -0.05, 0.05, 0.05, 0.05, 0.5,  0.75, -0.3, 1.0, 0, 0, 0, 0.221034, 0.813376, 0.918868},
		{-0.05, -0.05, -0.05, 0.05, 0.05, 0.05, 0.2,  0.75, 0.2,  1.0, 0, 0, 0, 0.421034, 0.813376, 0.918868},
		{-0.05, -0.05, -0.05, 0.05, 0.05, 0.05, -0.2, 0.75, -0.3, 1.0, 0, 0, 0, 0.221034, 0.413376, 0.918868},
		{-0.05, -0.05, -0.05, 0.05, 0.05, 0.05, -0.5, 0.75, -0.3, 1.0, 0, 0, 0, 0.221034, 0.413376, 0.418868}
	};

	for (int i = 0; i < 5; i++) {
		auto object = scene.add_objects();
		object->set_type(names[i] == "bin" ? Object_Type_BIN : Object_Type_BOX);
		object->set_id(names[i]);
		auto vi = values[i];
		auto pos = object->mutable_pos();
		pos->set_x((vi[0] + vi[3]) / 2 + vi[6]);
		pos->set_z((vi[1] + vi[4]) / 2 + vi[7]); // z and y are swapped here
		pos->set_y((vi[2] + vi[5]) / 2 + vi[8]); // because of the different coordinate systems
		auto size = object->mutable_size();
		size->set_width(vi[3] - vi[0]);
		size->set_height(vi[4] - vi[1]); // height and length are swapped here
		size->set_length(vi[5] - vi[2]); // because of the different coordinate systems

		object->mutable_orientation()->set_w(vi[9]);
		object->mutable_orientation()->set_x(vi[10]);
		object->mutable_orientation()->set_y(vi[11]);
		object->mutable_orientation()->set_z(vi[12]);

		object->mutable_color()->set_r(vi[13]);
		object->mutable_color()->set_g(vi[14]);
		object->mutable_color()->set_b(vi[15]);
	}
	return scene;
}

Selection vr_cobotics::obtainSelection(int box_id)
{
	Selection selection;
	selection.set_id(std::to_string(box_id));
	return selection;
}

void vr_cobotics::sender()
{
	nng::socket soc_pair_rec = nng::pair::open();
	soc_pair_rec.listen("tcp://127.0.0.1:6577");

	soc_pair_rec.dial("tcp://127.0.0.1:6576");

	Selection selection = obtainSelection(1);
	nng::view buf;
	int length = selection.ByteSize();
	void* data = nng_alloc(length);
	selection.SerializeToArray(data, length);

	buf = nng::view::view(data, length);
	soc_pair_rec.send(buf);
	nng::view rep_buf = soc_pair_rec.recv();

	// check the content
	if (rep_buf != "") {
		clear_movable_boxes();
		Scene scene;
		scene.ParseFromArray(rep_buf.data(), rep_buf.size());
		std::cout << "the number of objects: " << scene.objects_size() << std::endl;
		vec3 minp, maxp, trans;
		quat rot;
		rgb clr;
		for (auto& object : scene.objects())
		{
			std::cout << "type: " << object.type() << " name: " << object.id() << std::endl;
			if (object.type() == 1)
			{
				minp.x() = -object.size().length() / 2;
				minp.y() = -object.size().height() / 2;
				minp.z() = -object.size().width() / 2;
				maxp.x() = object.size().length() / 2;
				maxp.y() = object.size().height() / 2;
				maxp.z() = object.size().width() / 2;
				movable_boxes.emplace_back(minp, maxp);
				// exchange y and z, because ros uses a physical coordinate.
				trans.x() = object.pos().x();
				trans.y() = object.pos().z();
				trans.z() = object.pos().y();
				movable_box_translations.emplace_back(trans);
				rot.w() = object.orientation().w();
				rot.x() = object.orientation().x();
				rot.y() = object.orientation().z();
				rot.z() = object.orientation().y();
				movable_box_rotations.emplace_back(rot);
				clr.R() = object.color().r();
				clr.G() = object.color().g();
				clr.B() = object.color().b();
				movable_box_colors.emplace_back(clr);
				std::cout << object.pos().x() << std::endl;
			}
			else if (object.type() == 2)
			{
				is_trashbin = true;
				std::cout << "this is trash can" << std::endl;
			}
			else if (object.type() == 3)
			{
				std::cout << "this is robot arm" << std::endl;
			}
		}
		//save_boxes("boxes", movable_boxes, movable_box_colors, movable_box_translations, movable_box_rotations);
	}
}

void vr_cobotics::send_selection(int box_id)
{
	nng::socket soc_pair_rec = nng::pair::open();
	soc_pair_rec.listen("tcp://127.0.0.1:6577");

	soc_pair_rec.dial("tcp://127.0.0.1:6576");

	Selection selection = obtainSelection(box_id);
	nng::view buf;
	int length = selection.ByteSize();
	void* data = nng_alloc(length);
	selection.SerializeToArray(data, length);

	buf = nng::view::view(data, length);
	soc_pair_rec.send(buf);
	nng::view rep_buf = soc_pair_rec.recv();

	// check the content
	if (rep_buf != "") {
		Scene scene;
		scene.ParseFromArray(rep_buf.data(), rep_buf.size());
		std::cout << "the number of objects: " << scene.objects_size() << std::endl;
		vec3 minp, maxp, trans;
		quat rot;
		rgb clr;
		for (auto& object : scene.objects())
		{
			std::cout << "type: " << object.type() << " name: " << object.id() << std::endl;
			if (object.type() == 1)
			{
				minp.x() = -object.size().length() / 2;
				minp.y() = -object.size().height() / 2;
				minp.z() = -object.size().width() / 2;
				maxp.x() = object.size().length() / 2;
				maxp.y() = object.size().height() / 2;
				maxp.z() = object.size().width() / 2;
				movable_boxes.emplace_back(minp, maxp);
				// exchange y and z, because ros uses a physical coordinate.
				trans.x() = object.pos().x();
				trans.y() = object.pos().z();
				trans.z() = object.pos().y();
				movable_box_translations.emplace_back(trans);
				rot.w() = object.orientation().w();
				rot.x() = object.orientation().x();
				rot.y() = object.orientation().z();
				rot.z() = object.orientation().y();
				movable_box_rotations.emplace_back(rot);
				clr.R() = object.color().r();
				clr.G() = object.color().g();
				clr.B() = object.color().b();
				movable_box_colors.emplace_back(clr);
				std::cout << object.pos().x() << std::endl;
			}
		}
		//save_boxes("boxes", movable_boxes, movable_box_colors, movable_box_translations, movable_box_rotations);
	}
}

void vr_cobotics::receiver()
{

}

#include <cgv/base/register.h>

cgv::base::object_registration<vr_cobotics> vr_cobotics_reg("vr_cobotics");
