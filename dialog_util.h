#include <gtkmm.h>

bool valid_gps(Glib::ustring& text);
bool valid_number(Glib::ustring& text);

Glib::ustring gps_to_text(double latitude, double longitude);
Glib::ustring number_to_text(double value);
void text_to_gps(Glib::ustring& text, double &latitude, double &longitude);
void text_to_number(Glib::ustring& text, double &value);
