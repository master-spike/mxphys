#ifndef MXPHYS_EVENT_EVENTMANAGER_H
#define MXPHYS_EVENT_EVENTMANAGER_H

#include <unordered_map>
#include <vector>
#include <functional>
#include <memory>
#include <utility>

#include "event_types.h"
#include "event_data.h"

namespace mxphys::event {

class eventmanager {
private:
    typedef std::function<void(std::unique_ptr<event_data> const&)> event_listener;
    std::unordered_map<event_type, std::vector<event_listener>> m_registered_listeners;
    std::unordered_map<event_type, std::vector<std::unique_ptr<event_data>>> m_eventdata;
public:
    
    void add_event(event_type etype, std::unique_ptr<event_data>&& edata);
    void emit_all();
    void register_listener(event_type etype, const event_listener& listener);
};

}




#endif

