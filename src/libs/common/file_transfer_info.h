#include <string>
#include <cstdint>
#include <vector>
#include <memory>

class FileTransferInfo
{
public:
	//enum class TransferStatus :uint32_t { MANAGER_TO_SEND, MANAGER_TO_RETRIEVE, MANAGER_REQUESTED_OK, MANAGER_REQUESTED_FAILED, 
	//	MANAGER_RECEIVED_OK, MANAGER_RECEIVED_FAILED, MANAGER_SENT, WORKER_RECEIVED_OK, WORKER_RECEIVED_FAILED };
	//std::string filename_on_manager;
	enum class TransferType :uint32_t { SEND_TO_WORKERS, RETRIEVE_FROM_WORKER, UNKNWOWN };
	TransferType transfer_type;
	int file_number_on_manager;
	int file_number_on_worker;
	int target_run_id = -1;
	bool request_sent = false;
	bool file_received = false;
	std::vector<int> required_workers;
	//std::vector<int> workers_received_from;
	//std::vector<int> workers_sent_to;

private:

};
