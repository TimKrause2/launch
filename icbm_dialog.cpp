#include "icbm_dialog.h"
#include "dialog_util.h"
#include <ios>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <cmath>


#define RAD_PER_DEG (M_PI/180.0)
#define DEG_PER_RAD (180.0/M_PI)
#define ICBM_LAT_VEHICLE  48.935333180089664*RAD_PER_DEG
#define ICBM_LONG_VEHICLE -54.58524740913247*RAD_PER_DEG
#define ICBM_HGT_VEHICLE  5.0
#define ICBM_LAT_TARGET   35.94248310361987*RAD_PER_DEG
#define ICBM_LONG_TARGET  -5.617076200801353*RAD_PER_DEG
#define ICBM_HGT_TARGET   1000.0

ICBMLaunchData::ICBMLaunchData()
{
    lat_launch = ICBM_LAT_VEHICLE;
    long_launch = ICBM_LONG_VEHICLE;
    hgt_launch = ICBM_HGT_VEHICLE;
    lat_target = ICBM_LAT_TARGET;
    long_target = ICBM_LONG_TARGET;
    hgt_target = ICBM_HGT_TARGET;
}

ICBMDialog::ICBMDialog(Gtk::Window& parent, ICBMLaunchData &icbm_data)
    :
      icbm_data(icbm_data),
      m_LabelLaunch("Launch site coordinates (latitude, longitude)"),
      m_LabelLaunchHgt("Launch site height(m)"),
      m_LabelTarget("Target coordinates (latitude, longitude)"),
      m_LabelTargetHgt("Taraget height(m)"),
      m_ButtonBox(Gtk::Orientation::HORIZONTAL, 0),
      m_Button_OK("Schedule"),
      m_Button_Cancel("Cancel")
{
    set_transient_for(parent);
    set_destroy_with_parent(true);

    auto provider = new_entry_provider();

    entry_set_font(m_EntryLaunch, provider);
    entry_set_font(m_EntryLaunchHgt, provider);
    entry_set_font(m_EntryTarget, provider);
    entry_set_font(m_EntryTargetHgt, provider);

    set_title("Ballistic Launch");
    set_child(m_Grid);
    m_Grid.set_expand(true);

    m_Grid.attach(m_LabelLaunch, 0, 0);
    m_Grid.attach(m_EntryLaunch, 1, 0);
    m_EntryLaunch.set_text(gps_to_text(icbm_data.lat_launch, icbm_data.long_launch));
    m_EntryLaunch.set_hexpand(true);

    m_Grid.attach(m_LabelLaunchHgt, 0, 1);
    m_Grid.attach(m_EntryLaunchHgt, 1, 1);
    m_EntryLaunchHgt.set_text(number_to_text(icbm_data.hgt_launch));
    m_EntryLaunchHgt.set_hexpand(true);

    m_Grid.attach(m_LabelTarget, 0, 2);
    m_Grid.attach(m_EntryTarget, 1, 2);
    m_EntryTarget.set_text(gps_to_text(icbm_data.lat_target, icbm_data.long_target));
    m_EntryTarget.set_hexpand(true);

    m_Grid.attach(m_LabelTargetHgt, 0, 3);
    m_Grid.attach(m_EntryTargetHgt, 1, 3);
    m_EntryTargetHgt.set_text(number_to_text(icbm_data.hgt_target));
    m_EntryTargetHgt.set_hexpand(true);

    m_Grid.attach(m_ButtonBox, 0, 4, 2, 1);
    m_ButtonBox.set_halign(Gtk::Align::END);
    m_ButtonBox.append(m_Button_OK);
    m_ButtonBox.append(m_Button_Cancel);

    m_EntryLaunch.signal_activate().connect(sigc::mem_fun(*this,&ICBMDialog::on_launch_activate));
    m_EntryLaunchHgt.signal_activate().connect(sigc::mem_fun(*this,&ICBMDialog::on_launch_hgt_activate));
    m_EntryTarget.signal_activate().connect(sigc::mem_fun(*this,&ICBMDialog::on_target_activate));
    m_EntryTargetHgt.signal_activate().connect(sigc::mem_fun(*this,&ICBMDialog::on_target_hgt_activate));
}

ICBMDialog::~ICBMDialog()
{
}

void ICBMDialog::buttons_clicked_connect(
        const sigc::slot<void(const Glib::ustring&)>& slot)
{
    m_Button_OK.signal_clicked().connect(sigc::bind(slot, "OK"));
    m_Button_Cancel.signal_clicked().connect(sigc::bind(slot, "Cancel"));
}

void ICBMDialog::on_launch_activate()
{
    auto text = m_EntryLaunch.get_text();
    if(valid_gps(text)){
        text_to_gps(text, icbm_data.lat_launch, icbm_data.long_launch);
    }else{
        gpsAlertDialog(*this, "Launch site");
        m_EntryLaunch.grab_focus();
    }
}

void ICBMDialog::on_launch_hgt_activate()
{
    auto text = m_EntryLaunchHgt.get_text();
    if(valid_number(text)){
        text_to_number(text, icbm_data.hgt_launch);
    }else{
        numberAlertDialog(*this, "Launch site height.");
        m_EntryLaunchHgt.grab_focus();
    }
}

void ICBMDialog::on_target_activate()
{
    auto text = m_EntryTarget.get_text();
    if(valid_gps(text)){
        text_to_gps(text, icbm_data.lat_target, icbm_data.long_target);
    }else{
        gpsAlertDialog(*this, "Target");
        m_EntryTarget.grab_focus();
    }
}

void ICBMDialog::on_target_hgt_activate()
{
    auto text = m_EntryTargetHgt.get_text();
    if(valid_number(text)){
        text_to_number(text, icbm_data.hgt_target);
    }else{
        numberAlertDialog(*this, "Target height.");
        m_EntryTargetHgt.grab_focus();
    }
}
